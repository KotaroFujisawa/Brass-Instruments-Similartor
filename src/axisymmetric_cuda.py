"""CUDA Driver + NVRTC host code. The timestep graph remains on the device."""
import ctypes as ct
import os
from pathlib import Path
import sys

from axisymmetric_mesh import loss_fit
from trumpet import C, ROOT, np

PTR, DEV, INT, UINT, SIZE = ct.c_void_p, ct.c_uint64, ct.c_int, ct.c_uint, ct.c_size_t


class Config(ct.Structure):
    _fields_ = [(key,INT) for key in ('nc','nf','inlet','bell','stride','mode','feedback',
        'thermal_poles','radiation_order','reservoir_inflow')]+[
        (key,ct.c_float) for key in ('dt','delta','area_in','area_bell','r_bell','guard','viscous_dc')]+[
        ('filters_b',ct.c_float*5),('filters_w',ct.c_float*5)]+[
        (key,ct.c_double) for key in ('lip_hz','force','rest_y','hold','attack','release',
            'lip_width','lip_thickness','lip_joint','lip_rest_x','lip_mass_product','lip_cup_area','lip_q',
            'attack_boost','attack_decay','lip_gap_frequency_ratio','lip_gap_cubic_per_m2')]


TRUMPET_LIPS = dict(width=.007, thickness=.002, joint=.004, rest_x=.001,
                    mass_product=1.5, cup_area=.00023, quality=3.)


def load_nvrtc():
    if sys.platform == 'win32':
        directories = [ROOT/'tmp/cuda_runtime/nvidia/cuda_nvrtc/bin']
        directories += [Path(p)/'nvidia/cuda_nvrtc/bin' for p in sys.path if p]
        if os.environ.get('CUDA_NVRTC_DIR'):
            directories.insert(0,Path(os.environ['CUDA_NVRTC_DIR']))
        if os.environ.get('CUDA_PATH'):
            directories.append(Path(os.environ['CUDA_PATH'])/'bin')
        for directory in directories:
            dlls = sorted(directory.glob('nvrtc64_*.dll'))
            if dlls:
                handle = os.add_dll_directory(str(directory.resolve()))
                builtins = [ct.WinDLL(str(p.resolve())) for p in directory.glob('nvrtc-builtins*.dll')]
                library = ct.CDLL(str(dlls[-1].resolve()))
                library.dll_directory = handle
                library.builtins = builtins
                return library
        raise RuntimeError('NVRTC DLL not found. Install nvidia-cuda-nvrtc-cu12 or set CUDA_NVRTC_DIR.')
    for name in ('libnvrtc.so.12','libnvrtc.so'):
        try:
            return ct.CDLL(name)
        except OSError:
            continue
    raise RuntimeError('Set LD_LIBRARY_PATH to the CUDA NVRTC library directory.')


def compile_cuda(source, major, minor):
    lib=load_nvrtc()
    prog=PTR()
    lib.nvrtcCreateProgram.argtypes=[ct.POINTER(PTR),ct.c_char_p,ct.c_char_p,INT,PTR,PTR]
    lib.nvrtcCompileProgram.argtypes=[PTR,INT,ct.POINTER(ct.c_char_p)]
    for name in ('nvrtcGetProgramLogSize','nvrtcGetCUBINSize'):
        getattr(lib,name).argtypes=[PTR,ct.POINTER(SIZE)]
    for name in ('nvrtcGetProgramLog','nvrtcGetCUBIN'):
        getattr(lib,name).argtypes=[PTR,PTR]
    lib.nvrtcDestroyProgram.argtypes=[ct.POINTER(PTR)]
    if lib.nvrtcCreateProgram(ct.byref(prog),source,b'axisymmetric_euler.cu',0,None,None):
        raise RuntimeError('nvrtcCreateProgram failed')
    try:
        options=(ct.c_char_p*3)(b'--std=c++17',f'--gpu-architecture=sm_{major}{minor}'.encode(),b'--restrict')
        code=lib.nvrtcCompileProgram(prog,len(options),options)
        if code:
            size=SIZE();lib.nvrtcGetProgramLogSize(prog,ct.byref(size))
            log=ct.create_string_buffer(size.value);lib.nvrtcGetProgramLog(prog,log)
            raise RuntimeError(log.value.decode(errors='replace'))
        size=SIZE()
        if lib.nvrtcGetCUBINSize(prog,ct.byref(size)):
            raise RuntimeError('nvrtcGetCUBINSize failed')
        result=ct.create_string_buffer(size.value)
        if lib.nvrtcGetCUBIN(prog,result):
            raise RuntimeError('nvrtcGetCUBIN failed')
        return result
    finally:
        lib.nvrtcDestroyProgram(ct.byref(prog))


class Device:
    def __init__(self):
        self.lib=ct.WinDLL('nvcuda.dll') if sys.platform=='win32' else ct.CDLL('libcuda.so.1')
        signatures={
            'cuInit':[UINT], 'cuDeviceGet':[ct.POINTER(INT),INT],
            'cuDeviceGetAttribute':[ct.POINTER(INT),INT,INT],
            'cuDeviceGetName':[PTR,INT,INT], 'cuDevicePrimaryCtxRetain':[ct.POINTER(PTR),INT],
            'cuCtxSetCurrent':[PTR], 'cuModuleLoadData':[ct.POINTER(PTR),PTR],
            'cuModuleGetFunction':[ct.POINTER(PTR),PTR,ct.c_char_p],
            'cuMemAlloc_v2':[ct.POINTER(DEV),SIZE], 'cuMemFree_v2':[DEV],
            'cuMemcpyHtoD_v2':[DEV,PTR,SIZE], 'cuMemcpyDtoH_v2':[PTR,DEV,SIZE],
            'cuStreamCreate':[ct.POINTER(PTR),UINT], 'cuStreamSynchronize':[PTR],
            'cuStreamDestroy_v2':[PTR], 'cuStreamBeginCapture_v2':[PTR,INT],
            'cuStreamEndCapture':[PTR,ct.POINTER(PTR)],
            'cuGraphInstantiateWithFlags':[ct.POINTER(PTR),PTR,ct.c_uint64],
            'cuGraphDestroy':[PTR], 'cuGraphExecDestroy':[PTR], 'cuGraphLaunch':[PTR,PTR],
            'cuLaunchKernel':[PTR,UINT,UINT,UINT,UINT,UINT,UINT,UINT,PTR,PTR,PTR],
            'cuGetErrorName':[INT,ct.POINTER(ct.c_char_p)]}
        for name,signature in signatures.items():
            fn=getattr(self.lib,name);fn.argtypes=signature;fn.restype=INT
        self.call('cuInit',0)
        dev=INT();self.call('cuDeviceGet',ct.byref(dev),0)
        label=ct.create_string_buffer(128);self.call('cuDeviceGetName',label,128,dev)
        self.name=label.value.decode()
        self.context=PTR();self.call('cuDevicePrimaryCtxRetain',ct.byref(self.context),dev)
        self.call('cuCtxSetCurrent',self.context)
        capability=[]
        for attr in [75,76]:
            value=INT();self.call('cuDeviceGetAttribute',ct.byref(value),attr,dev)
            capability.append(value.value)
        blob=compile_cuda((ROOT/'src/axisymmetric_euler.cu').read_bytes(),*capability)
        self.module=PTR();self.call('cuModuleLoadData',ct.byref(self.module),blob)
        self.flux=self.kernel('fluxes');self.advance=self.kernel('advance')

    def call(self,name,*args):
        result=getattr(self.lib,name)(*args)
        if result:
            error=ct.c_char_p();self.lib.cuGetErrorName(result,ct.byref(error))
            raise RuntimeError(f'{name}: {error.value} ({result})')

    def kernel(self,name):
        fn=PTR();self.call('cuModuleGetFunction',ct.byref(fn),self.module,name.encode())
        return fn


class Run:
    def __init__(self, device, mesh, seconds, *, lip_hz=490., pressure=7600., rest_y=0.,
                 hold=None, attack=.045, release=.065, mode=0, feedback=True, lip_parameters=None,
                 attack_boost=0., attack_decay=.024, thermoviscous=False, radiation_order=1,
                 reservoir_inflow=False,
                 initial=None, trace_fs=192000, cfl=.73, guard=1.58, chunk=256):
        self.device=device;self.mesh=mesh;self.trace_fs=trace_fs
        self.chunk=chunk
        mesh_cells,faces=mesh['cells'],mesh['faces']
        nc,nf=len(mesh_cells),len(faces)
        if max(mesh['rings'][[0,-1]])>128:
            raise ValueError('The endpoint meshes may contain at most 128 rings.')
        self.stride=int(np.ceil(C*guard/(mesh['length_scale']*cfl*trace_fs)))
        self.dt=1/(trace_fs*self.stride)
        self.records=round(seconds*trace_fs);self.steps=self.records*self.stride
        self.progress=0
        if hold is None:hold=seconds-.14
        if attack<=0 or attack_boost<0 or attack_decay<=0:
            raise ValueError('Attack time and decay must be positive, and the pressure boost nonnegative.')
        lips = TRUMPET_LIPS if lip_parameters is None else TRUMPET_LIPS | lip_parameters
        if radiation_order not in (1,2):
            raise ValueError('Radiation filter order must be 1 or 2.')
        if thermoviscous:
            from thermoviscous import KINEMATIC_VISCOSITY, memory_coefficients
            loss_coefficients,self.loss_error=memory_coefficients(mesh_cells['radius2'].astype('f8'),self.dt)
            poles=len(loss_coefficients)
            viscous_dc=8*KINEMATIC_VISCOSITY*self.dt
            theta=weights=np.zeros(5)
        else:
            loss_coefficients=np.zeros((1,4))
            poles=0
            viscous_dc=0.
            theta,weights,self.loss_error=loss_fit()
        config=Config(nc,nf,mesh['inlet'],mesh['bell'],self.stride,mode,feedback,
                      poles,radiation_order,reservoir_inflow,
                      self.dt,C*self.dt,np.pi*mesh['radius'][0]**2,np.pi*mesh['radius'][-1]**2,
                      mesh['radius'][-1],guard,viscous_dc,(ct.c_float*5)(*(-np.expm1(-theta*self.dt))),
                      (ct.c_float*5)(*(weights*self.dt)),lip_hz,pressure,rest_y,hold,attack,release,
                      *(lips[key] for key in TRUMPET_LIPS),attack_boost,attack_decay,
                      lips.get('gap_frequency_ratio',1.),lips.get('gap_cubic_per_m2',0.))
        self.buffers=[]
        self.stream=PTR();device.call('cuStreamCreate',ct.byref(self.stream),1)
        self.graph=PTR();self.executable=PTR()
        nc4=(nc,4)
        self.state=self.upload(np.zeros(nc4,dtype='f4') if initial is None else np.asarray(initial,dtype='f4'))
        self.predicted=self.upload(np.zeros(nc4,dtype='f4'))
        self.cells=self.upload(mesh_cells);self.faces=self.upload(faces)
        self.flux=self.upload(np.zeros((nf,4),dtype='f4'))
        self.memory=self.upload(np.zeros((poles,nc,4) if poles else (5,nc,2),dtype='f4'))
        self.loss_coefficients=self.upload(loss_coefficients.astype('f4'))
        self.config=self.upload(np.frombuffer(bytes(config),dtype='u1'))
        self.clock=self.upload(np.zeros(4,dtype='i4'))
        boundary=np.zeros(15,dtype='f8');boundary[:4]=[lips['rest_x'],rest_y,lips['rest_x'],rest_y]
        self.boundary=self.upload(boundary)
        self.trace=self.upload(np.zeros((self.records,11),dtype='f4'))
        self.kernels=[]
        for stage in (0,1):
            self.kernels.append((device.flux,(mesh['inlet']+127)//128+2,self.args(
                self.state,self.predicted,self.cells,self.faces,self.flux,self.config,self.clock,self.boundary,INT(stage))))
            self.kernels.append((device.advance,(nc+127)//128,self.args(self.state,self.predicted,self.cells,
                self.flux,self.memory,self.loss_coefficients,self.config,self.clock,self.boundary,self.trace,INT(stage))))
        device.call('cuStreamBeginCapture_v2',self.stream,0)
        for _ in range(chunk):self.launch_step()
        device.call('cuStreamEndCapture',self.stream,ct.byref(self.graph))
        device.call('cuGraphInstantiateWithFlags',ct.byref(self.executable),self.graph,0)

    def upload(self,a):
        a=np.ascontiguousarray(a)
        p=DEV();self.device.call('cuMemAlloc_v2',ct.byref(p),a.nbytes)
        self.buffers.append(p)
        self.device.call('cuMemcpyHtoD_v2',p,a.ctypes.data,a.nbytes)
        return p

    def download(self,p,shape,dtype):
        a=np.empty(shape,dtype=dtype)
        self.device.call('cuMemcpyDtoH_v2',a.ctypes.data,p,a.nbytes)
        return a

    @staticmethod
    def args(*items):
        # Keep scalar storage live for the host API during stream capture.
        return ((PTR*len(items))(*(ct.addressof(v) for v in items)),items)

    def launch_step(self):
        for kernel,blocks,(args,_) in self.kernels:
            self.device.call('cuLaunchKernel',kernel,blocks,1,1,128,1,1,0,self.stream,args,None)

    def until(self, seconds=None):
        limit=self.steps if seconds is None else min(round(seconds*self.trace_fs)*self.stride,self.steps)
        if limit<self.progress:raise ValueError('Time cannot run backward')
        graphs,remainder=divmod(limit-self.progress,self.chunk)
        for i in range(graphs):
            self.device.call('cuGraphLaunch',self.executable,self.stream)
            if i%64==63:self.device.call('cuStreamSynchronize',self.stream)
        for _ in range(remainder):self.launch_step()
        self.device.call('cuStreamSynchronize',self.stream)
        self.progress=limit
        step,error,self.fallbacks,_=self.download(self.clock,(4,),'i4')
        if step!=limit or error:
            raise RuntimeError(f'Invalid solver state at step {step}/{limit}, status {error}; refine timestep or increase wave-speed reserve')

    def field(self):
        return self.download(self.state,(len(self.mesh['cells']),4),'f4')

    def history(self):
        return self.download(self.trace,(self.records,11),'f4')[:self.progress//self.stride]

    def __enter__(self):return self

    def __exit__(self,*exc):
        for p in reversed(self.buffers):self.device.call('cuMemFree_v2',p)
        self.device.call('cuGraphExecDestroy',self.executable)
        self.device.call('cuGraphDestroy',self.graph)
        self.device.call('cuStreamDestroy_v2',self.stream)
