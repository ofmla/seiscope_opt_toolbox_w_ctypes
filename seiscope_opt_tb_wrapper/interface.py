import os
import ctypes
from ctypes import POINTER, c_int, c_float, c_bool


class UserDefined(ctypes.Structure):
    """Creates a struct to match the optim Fortran derived type in
    SEISCOPE optimization toolbox. Fields of the derived type are
    stored in the _fields_ attribute, which is a dict.
    """
    _fields_ = [
        ('debug', c_bool),
        ('threshold', c_float),
        ('print_flag', c_int),
        ('first_ls', c_bool),
        ('task', c_int),
        ('nls_max', c_int),
        ('cpt_ls', c_int),
        ('nfwd_pb', c_int),
        ('cpt_iter', c_int),
        ('niter_max', c_int),
        ('f0', c_float),
        ('fk', c_float),
        ('conv', c_float),
        ('m1', c_float),
        ('m2', c_float),
        ('mult_factor', c_float),
        ('alpha_L', c_float),
        ('alpha_R', c_float),
        ('alpha', c_float),
        ('q0', c_float),
        ('q', c_float),
        ('cpt_lbfgs', c_int),
        ('l', c_int),
        ('conv_CG', c_bool),
        ('cpt_iter_CG', c_int),
        ('niter_max_CG', c_int),
        ('nhess', c_int),
        ('CG_phase', c_int),
        ('comm', c_int),
        ('qk_CG', c_float),
        ('qkm1_CG', c_float),
        ('hessian_term', c_float),
        ('eta', c_float),
        ('norm_grad', c_float),
        ('norm_grad_m1', c_float),
        ('norm_residual', c_float),
        ('dHd', c_float),
        ('res_scal_respreco', c_float),
        ('alpha_CG', c_float)
    ]

    def __repr__(self):
        """Print a representation of the derived type."""
        template = (
            'UserDefined(debug={self.debug}, '
            'threshold={self.threshold}, '
            'print_flag={self.print_flag}, '
            'first_ls={self.first_ls}, '
            'task={self.task}, '
            'nls_max={self.nls_max}, '
            'cpt_ls={self.cpt_ls}, '
            'nfwd_pb={self.nfwd_pb}, '
            'cpt_iter={self.cpt_iter}, '
            'niter_max={self.niter_max}, '
            'f0={self.f0}, '
            'fk={self.fk}, '
            'conv={self.conv}, '
            'm1={self.m1}, '
            'm2={self.m2}, '
            'mult_factor={self.mult_factor}, '
            'alpha_L={self.alpha_L}, '
            'alpha_R={self.alpha_R}, '
            'alpha={self.alpha}, '
            'q0={self.q0}, '
            'q={self.q}, '
            'cpt_lbfgs={self.cpt_lbfgs}, '
            'l={self.l},'
            'conv_CG={self.conv_CG},'
            'cpt_iter_CG={self.cpt_iter_CG},'
            'niter_max_CG={self.niter_max_CG},'
            'nhess={self.nhess},'
            'CG_phase={self.CG_phase},'
            'comm={self.comm},'
            'qk_CG={self.qk_CG},'
            'qkm1_CG={self.qkm1_CG},'
            'hessian_term={self.hessian_term},'
            'eta={self.eta},'
            'norm_grad={self.norm_grad},'
            'norm_grad_m1_CG={self.norm_grad_m1},'
            'norm_residual={self.norm_residual},'
            'dHd={self.dHd},'
            'res_scal_respreco={self.res_scal_respreco},'
            'alpha_CG={self.alpha_CG})'
        )
        return template.format(self=self)


# Get the location of the shared library file.
here = os.path.dirname(os.path.abspath(__file__))
lib_file = os.path.join(here, '..', 'lib', 'libOPTIM.so')

# This is how a dll/so library is loaded
lib_sotb = ctypes.cdll.LoadLibrary(lib_file)


class sotb_wrapper(object):

    _ctypes_pstd = lib_sotb['PSTD']
    _ctypes_pnlcg = lib_sotb['PNLCG']
    _ctypes_plbfgs = lib_sotb['PLBFGS']
    _ctypes_lbfgs = lib_sotb['LBFGS']
    _ctypes_trn = lib_sotb['TRN']
    _ctypes_ptrn = lib_sotb['PTRN']

    _ctypes_pstd.argtypes = [POINTER(c_int), POINTER(c_float), POINTER(c_float),
                             POINTER(c_float), POINTER(c_float), POINTER(UserDefined),
                             POINTER(c_int), POINTER(c_float), POINTER(c_float)]
    _ctypes_pstd.restype = None

    _ctypes_pnlcg.argtypes = [POINTER(c_int), POINTER(c_float), POINTER(c_float),
                              POINTER(c_float), POINTER(c_float), POINTER(UserDefined),
                              POINTER(c_int), POINTER(c_float), POINTER(c_float)]
    _ctypes_pnlcg.restype = None

    _ctypes_plbfgs.argtypes = [POINTER(c_int), POINTER(c_float), POINTER(c_float),
                               POINTER(c_float), POINTER(c_float),
                               POINTER(UserDefined), POINTER(c_int), POINTER(c_float),
                               POINTER(c_float)]
    _ctypes_plbfgs.restype = None

    _ctypes_lbfgs.argtypes = [POINTER(c_int), POINTER(c_float), POINTER(c_float),
                              POINTER(c_float), POINTER(UserDefined), POINTER(c_int),
                              POINTER(c_float), POINTER(c_float)]
    _ctypes_lbfgs.restype = None

    _ctypes_trn.argtypes = [POINTER(c_int), POINTER(c_float), POINTER(c_float),
                            POINTER(c_float), POINTER(c_float), POINTER(c_float),
                            POINTER(UserDefined), POINTER(c_int),
                            POINTER(c_float), POINTER(c_float)]
    _ctypes_trn.restype = None

    _ctypes_ptrn.argtypes = [POINTER(c_int), POINTER(c_float), POINTER(c_float),
                             POINTER(c_float), POINTER(c_float), POINTER(c_float),
                             POINTER(c_float), POINTER(c_float), POINTER(UserDefined),
                             POINTER(c_int), POINTER(c_float), POINTER(c_float)]
    _ctypes_ptrn.restype = None

    def __init__(self):
        self.udf = UserDefined()

    # wrapping functions for each method
    def PSTD(self, n, x, fcost, grad, grad_preco, flag, lb=None, ub=None):
        lb_pass = lb.ctypes.data_as(POINTER(c_float)) if lb else lb
        ub_pass = ub.ctypes.data_as(POINTER(c_float)) if ub else ub
        self._ctypes_pstd(ctypes.byref(n), x.ctypes.data_as(POINTER(c_float)),
                          ctypes.byref(fcost), grad.ctypes.data_as(POINTER(c_float)),
                          grad_preco.ctypes.data_as(POINTER(c_float)),
                          ctypes.byref(self.udf), ctypes.byref(flag), lb_pass, ub_pass)

    def PNLCG(self, n, x, fcost, grad, grad_preco, flag, lb=None, ub=None):
        lb_pass = lb.ctypes.data_as(POINTER(c_float)) if lb else lb
        ub_pass = ub.ctypes.data_as(POINTER(c_float)) if ub else ub
        self._ctypes_pnlcg(ctypes.byref(n), x.ctypes.data_as(POINTER(c_float)),
                           ctypes.byref(fcost), grad.ctypes.data_as(POINTER(c_float)),
                           grad_preco.ctypes.data_as(POINTER(c_float)),
                           ctypes.byref(self.udf), ctypes.byref(flag), lb_pass, ub_pass)

    def PLBFGS(self, n, x, fcost, grad, grad_preco, flag, lb=None, ub=None):
        lb_pass = lb.ctypes.data_as(POINTER(c_float)) if lb else lb
        ub_pass = ub.ctypes.data_as(POINTER(c_float)) if ub else ub
        self._ctypes_plbfgs(ctypes.byref(n), x.ctypes.data_as(POINTER(c_float)),
                            ctypes.byref(fcost), grad.ctypes.data_as(POINTER(c_float)),
                            grad_preco.ctypes.data_as(POINTER(c_float)),
                            ctypes.byref(self.udf), ctypes.byref(flag), lb_pass, ub_pass)

    def LBFGS(self, n, x, fcost, grad, flag, lb=None, ub=None):
        lb_pass = lb.ctypes.data_as(POINTER(c_float)) if lb is not None else lb
        ub_pass = ub.ctypes.data_as(POINTER(c_float)) if ub is not None else ub
        self._ctypes_lbfgs(ctypes.byref(n), x.ctypes.data_as(POINTER(c_float)),
                           ctypes.byref(fcost), grad.ctypes.data_as(POINTER(c_float)),
                           ctypes.byref(self.udf), ctypes.byref(flag), lb_pass, ub_pass)

    def TRN(self, n, x, fcost, grad, d, hd, flag, lb=None, ub=None):
        lb_pass = lb.ctypes.data_as(POINTER(c_float)) if lb is not None else lb
        ub_pass = ub.ctypes.data_as(POINTER(c_float)) if ub is not None else ub
        self._ctypes_trn(ctypes.byref(n), x.ctypes.data_as(POINTER(c_float)),
                         ctypes.byref(fcost), grad.ctypes.data_as(POINTER(c_float)),
                         d.ctypes.data_as(POINTER(c_float)),
                         hd.ctypes.data_as(POINTER(c_float)), ctypes.byref(self.udf),
                         ctypes.byref(flag), lb_pass, ub_pass)

    def PTRN(self, n, x, fcost, grad, grad_preco, residual, residual_preco, d, hd, flag,
             lb=None, ub=None):
        lb_pass = lb.ctypes.data_as(POINTER(c_float)) if lb is not None else lb
        ub_pass = ub.ctypes.data_as(POINTER(c_float)) if ub is not None else ub
        self._ctypes_ptrn(ctypes.byref(n), x.ctypes.data_as(POINTER(c_float)),
                          ctypes.byref(fcost), grad.ctypes.data_as(POINTER(c_float)),
                          grad_preco.ctypes.data_as(POINTER(c_float)),
                          residual.ctypes.data_as(POINTER(c_float)),
                          residual_preco.ctypes.data_as(POINTER(c_float)),
                          d.ctypes.data_as(POINTER(c_float)),
                          hd.ctypes.data_as(POINTER(c_float)), ctypes.byref(self.udf),
                          ctypes.byref(flag), lb_pass, ub_pass)
