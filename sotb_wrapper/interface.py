"""Optimization methods code."""
from __future__ import annotations

import ctypes
from ctypes import byref, c_float, c_int, POINTER
from pathlib import Path
import sys
from typing import Optional

import numpy as np
from numpy.typing import ArrayLike, NDArray


class UserDefined(ctypes.Structure):
    """Creates python object representing a C struct.

    Creates a struct to match the optim Fortran derived type in SEISCOPE optimization
    toolbox. Fields of the derived type are stored in the _fields_ attribute, which is
    a dict.
    """

    _fields_ = [
        ("debug", c_int),
        ("threshold", c_float),
        ("print_flag", c_int),
        ("first_ls", c_int),
        ("task", c_int),
        ("nls_max", c_int),
        ("cpt_ls", c_int),
        ("nfwd_pb", c_int),
        ("cpt_iter", c_int),
        ("niter_max", c_int),
        ("f0", c_float),
        ("fk", c_float),
        ("conv", c_float),
        ("m1", c_float),
        ("m2", c_float),
        ("mult_factor", c_float),
        ("alpha_L", c_float),
        ("alpha_R", c_float),
        ("alpha", c_float),
        ("q0", c_float),
        ("q", c_float),
        ("cpt_lbfgs", c_int),
        ("l", c_int),
        ("conv_CG", c_int),
        ("cpt_iter_CG", c_int),
        ("niter_max_CG", c_int),
        ("nhess", c_int),
        ("CG_phase", c_int),
        ("comm", c_int),
        ("qk_CG", c_float),
        ("qkm1_CG", c_float),
        ("hessian_term", c_float),
        ("eta", c_float),
        ("norm_grad", c_float),
        ("norm_grad_m1", c_float),
        ("norm_residual", c_float),
        ("dHd", c_float),
        ("res_scal_respreco", c_float),
        ("alpha_CG", c_float),
    ]

    def __repr__(self) -> str:
        """Print a representation of the C struc (Fortran derived type).

        Returns:
            template: A string representation of an instance of `UserDefined` class
        """
        template = (
            "UserDefined(debug={self.debug}, "
            "threshold={self.threshold}, "
            "print_flag={self.print_flag}, "
            "first_ls={self.first_ls}, "
            "task={self.task}, "
            "nls_max={self.nls_max}, "
            "cpt_ls={self.cpt_ls}, "
            "nfwd_pb={self.nfwd_pb}, "
            "cpt_iter={self.cpt_iter}, "
            "niter_max={self.niter_max}, "
            "f0={self.f0}, "
            "fk={self.fk}, "
            "conv={self.conv}, "
            "m1={self.m1}, "
            "m2={self.m2}, "
            "mult_factor={self.mult_factor}, "
            "alpha_L={self.alpha_L}, "
            "alpha_R={self.alpha_R}, "
            "alpha={self.alpha}, "
            "q0={self.q0}, "
            "q={self.q}, "
            "cpt_lbfgs={self.cpt_lbfgs}, "
            "l={self.l}, "
            "conv_CG={self.conv_CG}, "
            "cpt_iter_CG={self.cpt_iter_CG}, "
            "niter_max_CG={self.niter_max_CG}, "
            "nhess={self.nhess}, "
            "CG_phase={self.CG_phase}, "
            "comm={self.comm}, "
            "qk_CG={self.qk_CG}, "
            "qkm1_CG={self.qkm1_CG}, "
            "hessian_term={self.hessian_term}, "
            "eta={self.eta}, "
            "norm_grad={self.norm_grad}, "
            "norm_grad_m1_CG={self.norm_grad_m1}, "
            "norm_residual={self.norm_residual}, "
            "dHd={self.dHd}, "
            "res_scal_respreco={self.res_scal_respreco}, "
            "alpha_CG={self.alpha_CG})"
        )
        return template.format(self=self)


# Load shared library
path = Path(__file__).parent
name = "libsotb"
if sys.platform.startswith("linux"):
    name = f"{name}.so"
elif sys.platform == "darwin":
    name = f"{name}.dylib"
else:
    raise ImportError(f"Your OS is not supported: {sys.platform}")

# This is how a dll/so library is loaded
try:
    lib_sotb = ctypes.CDLL(str(path / name), winmode=0)

except Exception as e:
    import logging

    logger = logging.Logger("catch_all")
    logger.error(e, exc_info=True)
    github_url = "https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes#compiling"
    print()
    print(
        "Failed to load the required SEISCOPE OPTIMIZATION TOOLBOX (sotb) shared"
        + " library."
        + "\n"
        + "You are likely to see this message as consequence of the attempt to "
        + "run the tests "
        + "\n"
        + "from the local project directory (the one obtained with the clone command). "
        + "\n"
        + "\n"
        + "If your objective is to contribute to sotb-wrapper and you got a development "
        + "\n"
        + "copy via `git clone` command, you must have the nox package to be able to run"
        + "\n"
        + "the tests from the local repository. So an immediate step after cloning the "
        + "repo"
        + "\n"
        + "is to install nox in your python environment"
        + "\n"
        + "\n"
        + "python3 -m pip install nox"
        + "\n"
        + "\n"
        + "You can then make changes and use the command `nox -s tests` to execute the "
        + "tests session. "
        + "\n"
    )
    print()


class sotb_wrapper(object):
    """Includes all optimization methods from SEISCOPE Optimization toolbox."""

    _ctypes_pstd = lib_sotb["PSTD"]
    _ctypes_pnlcg = lib_sotb["PNLCG"]
    _ctypes_plbfgs = lib_sotb["PLBFGS"]
    _ctypes_lbfgs = lib_sotb["LBFGS"]
    _ctypes_trn = lib_sotb["TRN"]
    _ctypes_ptrn = lib_sotb["PTRN"]
    _ctypes_init = lib_sotb["set_inputs"]
    _ctypes_rosenbrock = lib_sotb["rosenbrock"]
    _ctypes_rosenbrock_hess = lib_sotb["rosenbrock_hess"]

    _ctypes_rosenbrock.argtypes = [
        POINTER(c_float),
        POINTER(c_float),
        POINTER(c_float),
    ]
    _ctypes_rosenbrock.restype = None

    _ctypes_rosenbrock_hess.argtypes = [
        POINTER(c_float),
        POINTER(c_float),
        POINTER(c_float),
    ]
    _ctypes_rosenbrock_hess.restype = None

    _ctypes_init.argtypes = [
        POINTER(UserDefined),
        c_float,
        c_int,
        POINTER(c_float),
        POINTER(c_int),
        POINTER(c_float),
        POINTER(c_int),
        POINTER(c_int),
        POINTER(c_int),
        POINTER(c_int),
    ]
    _ctypes_init.restype = None

    _ctypes_pstd.argtypes = [
        c_int,
        POINTER(c_float),
        POINTER(c_float),
        POINTER(c_float),
        POINTER(c_float),
        POINTER(UserDefined),
        POINTER(c_int),
        POINTER(c_float),
        POINTER(c_float),
    ]
    _ctypes_pstd.restype = None

    _ctypes_pnlcg.argtypes = [
        c_int,
        POINTER(c_float),
        POINTER(c_float),
        POINTER(c_float),
        POINTER(c_float),
        POINTER(UserDefined),
        POINTER(c_int),
        POINTER(c_float),
        POINTER(c_float),
    ]
    _ctypes_pnlcg.restype = None

    _ctypes_plbfgs.argtypes = [
        c_int,
        POINTER(c_float),
        POINTER(c_float),
        POINTER(c_float),
        POINTER(c_float),
        POINTER(c_float),
        POINTER(UserDefined),
        POINTER(c_int),
        POINTER(c_float),
        POINTER(c_float),
    ]
    _ctypes_plbfgs.restype = None

    _ctypes_lbfgs.argtypes = [
        c_int,
        POINTER(c_float),
        POINTER(c_float),
        POINTER(c_float),
        POINTER(UserDefined),
        POINTER(c_int),
        POINTER(c_float),
        POINTER(c_float),
    ]
    _ctypes_lbfgs.restype = None

    _ctypes_trn.argtypes = [
        c_int,
        POINTER(c_float),
        POINTER(c_float),
        POINTER(c_float),
        POINTER(c_float),
        POINTER(c_float),
        POINTER(UserDefined),
        POINTER(c_int),
        POINTER(c_float),
        POINTER(c_float),
    ]
    _ctypes_trn.restype = None

    _ctypes_ptrn.argtypes = [
        c_int,
        POINTER(c_float),
        POINTER(c_float),
        POINTER(c_float),
        POINTER(c_float),
        POINTER(c_float),
        POINTER(c_float),
        POINTER(c_float),
        POINTER(UserDefined),
        POINTER(c_int),
        POINTER(c_float),
        POINTER(c_float),
    ]
    _ctypes_ptrn.restype = None

    def __init__(self) -> None:
        """Instantiate a sotb_wrapper object with ``UserDefined`` property."""
        self.udf = UserDefined()

    # wrapping functions for each method
    def rosenbrock(self, x: ArrayLike) -> tuple[NDArray, float]:
        """Rosenbrock's function.

        Args:
            x: 1-D array of points at which the Rosenbrock function is to be computed

        Returns:
            grad: Gradient of the Rosenbrock function
            fcost: The value of the Rosenbrock function
        """
        fcost_ = c_float()
        grad = np.zeros(2, dtype=np.float32)
        self._ctypes_rosenbrock(
            x.ctypes.data_as(POINTER(c_float)),
            byref(fcost_),
            grad.ctypes.data_as(POINTER(c_float)),
        )
        return grad, float(fcost_.value)

    def rosenbrock_hess(self, x: ArrayLike, d: ArrayLike) -> NDArray:
        """Calculates the product of `d` and the Hessian of the Rosenbrock's function.

        Args:
            x: 1-D array of points at which the Hessian of the Rosenbrock
                function is to be computed
            d: Vector used in the conjugate gradient process by (P)TRN method

        Returns:
            hd: Product of `d` Hessian of the Rosenbrock function
        """
        hd = np.zeros(2, dtype=np.float32)
        self._ctypes_rosenbrock_hess(
            x.ctypes.data_as(POINTER(c_float)),
            d.ctypes.data_as(POINTER(c_float)),
            hd.ctypes.data_as(POINTER(c_float)),
        )
        return hd

    def set_inputs(
        self,
        fcost: float,
        niter_max: int,
        alpha: Optional[float] = None,
        nls_max: Optional[int] = None,
        conv: Optional[float] = None,
        print_flag: Optional[int] = None,
        l: Optional[int] = None,
        niter_max_CG: Optional[int] = None,
        debug: Optional[int] = None,
    ) -> None:
        """Initializes fields in UserDefined class.

        The function must be called after calculating the `fcost` for the initial guess
        and before the optimization method can be used.
        The following inputs are required: fcost, niter_max. For the others, if they
        are not present, then the default values are used

        Args:
            fcost: Cost function value associated with the initial guess
            niter_max: Maximum number of iterations
            alpha: The value to use as initial steplength. Defaults to 1.0
            nls_max: Maximum number of linesearch iterations. Defaults to 20.
            conv: Tolerance for the stopping criterion. By default, it is set to machine
                epsilon related to arithmetic with single precision real.
            print_flag: Specifies whether output files with info about the optimization
                process should be created. Accepted values are 1(default) and 0.
            l: Maximum number of stored pairs used for the l-BFGS approximation.
                Defaults to 10.
            niter_max_CG: Maximum number of inner conjugate gradient iterations.
                Defaults to 5.
            debug: Determines whether the highest level of details for output files needs
                to be considered. By default, debug is False, in which case the basic
                level of detail is used.
        """
        self._ctypes_init(
            byref(self.udf),
            c_float(fcost),
            c_int(niter_max),
            None if alpha is None else byref(c_float(alpha)),
            None if nls_max is None else byref(c_int(nls_max)),
            None if conv is None else byref(c_float(conv)),
            None if print_flag is None else byref(c_int(print_flag)),
            None if l is None else byref(c_int(l)),
            None if niter_max_CG is None else byref(c_int(niter_max_CG)),
            None if debug is None else byref(c_int(debug)),
        )

    def PSTD(
        self,
        n: int,
        x: ArrayLike,
        fcost: float,
        grad: ArrayLike,
        grad_preco: ArrayLike,
        flag: int,
        lb: Optional[ArrayLike] = None,
        ub: Optional[ArrayLike] = None,
    ) -> int:
        """Preconditioned steepest descent algorithm.

        Args:
            n: The parameter space dimension (size of the x vector)
            x: Solution vector
            fcost: Cost function value associated with x
            grad: Gradient at x
            grad_preco: Preconditioned gradient at x
            flag: Variable that indicates to the user on return what action he has to do,
                or the state of the algorithm
            lb: Vector of lower bounds for x. By default None
            ub: Vector of upper bounds for x. By default None

        Returns:
            flag: Variable that indicates to the user on return what action he has to do,
                or the state of the algorithm
        """
        lb_pass = lb.ctypes.data_as(POINTER(c_float)) if lb is not None else lb
        ub_pass = ub.ctypes.data_as(POINTER(c_float)) if ub is not None else ub
        cflag = c_int(flag)
        self._ctypes_pstd(
            c_int(n),
            x.ctypes.data_as(POINTER(c_float)),
            byref(c_float(fcost)),
            grad.ctypes.data_as(POINTER(c_float)),
            grad_preco.ctypes.data_as(POINTER(c_float)),
            byref(self.udf),
            byref(cflag),
            lb_pass,
            ub_pass,
        )
        return cflag.value

    def PNLCG(
        self,
        n: int,
        x: ArrayLike,
        fcost: float,
        grad: ArrayLike,
        grad_preco: ArrayLike,
        flag: int,
        lb: Optional[ArrayLike] = None,
        ub: Optional[ArrayLike] = None,
    ) -> int:
        """Preconditioned nonlinear conjugate gradient algorithm.

        Args:
            n: The parameter space dimension (size of the x vector)
            x: Solution vector
            fcost: Cost function value associated with x
            grad: Gradient at x
            grad_preco: Preconditioned gradient at x
            flag: Variable that indicates to the user on return what action he has to do,
                or the state of the algorithm
            lb: Vector of lower bounds for x. By default None
            ub: Vector of upper bounds for x. By default None

        Returns:
            flag: Variable that indicates to the user on return what action he has to do,
                or the state of the algorithm
        """
        lb_pass = lb.ctypes.data_as(POINTER(c_float)) if lb is not None else lb
        ub_pass = ub.ctypes.data_as(POINTER(c_float)) if ub is not None else ub
        cflag = c_int(flag)
        self._ctypes_pnlcg(
            n,
            x.ctypes.data_as(POINTER(c_float)),
            byref(c_float(fcost)),
            grad.ctypes.data_as(POINTER(c_float)),
            grad_preco.ctypes.data_as(POINTER(c_float)),
            byref(self.udf),
            byref(cflag),
            lb_pass,
            ub_pass,
        )
        return cflag.value

    def PLBFGS(
        self,
        n: int,
        x: ArrayLike,
        fcost: float,
        grad: ArrayLike,
        grad_preco: ArrayLike,
        q_plb: ArrayLike,
        flag: int,
        lb: Optional[ArrayLike] = None,
        ub: Optional[ArrayLike] = None,
    ) -> int:
        """Preconditioned l-BFGS algorithm.

        Args:
            n: The parameter space dimension (size of the x vector)
            x: Solution vector
            fcost: Cost function value associated with x
            grad: Gradient at x
            grad_preco: Preconditioned gradient at x
            q_plb: l-BFGS approximation of the inverse Hessian operator
            flag: Variable that indicates to the user on return what action he has to do,
                or the state of the algorithm
            lb: Vector of lower bounds for x. By default None
            ub: Vector of upper bounds for x. By default None

        Returns:
            flag: Variable that indicates to the user on return what action he has to do,
                or the state of the algorithm
        """
        lb_pass = lb.ctypes.data_as(POINTER(c_float)) if lb is not None else lb
        ub_pass = ub.ctypes.data_as(POINTER(c_float)) if ub is not None else ub
        cflag = c_int(flag)
        self._ctypes_plbfgs(
            n,
            x.ctypes.data_as(POINTER(c_float)),
            byref(c_float(fcost)),
            grad.ctypes.data_as(POINTER(c_float)),
            grad_preco.ctypes.data_as(POINTER(c_float)),
            q_plb.ctypes.data_as(POINTER(c_float)),
            byref(self.udf),
            byref(cflag),
            lb_pass,
            ub_pass,
        )
        return cflag.value

    def LBFGS(
        self,
        n: int,
        x: ArrayLike,
        fcost: float,
        grad: ArrayLike,
        flag: int,
        lb: Optional[ArrayLike] = None,
        ub: Optional[ArrayLike] = None,
    ) -> int:
        """l-BFGS algorithm.

        Args:
            n: The parameter space dimension (size of the x vector)
            x: Solution vector
            fcost: Cost function value associated with x
            grad: Gradient at x
            flag: Variable that indicates to the user on return what action he has to do,
                or the state of the algorithm
            lb: Vector of lower bounds for x. By default None
            ub: Vector of upper bounds for x. By default None

        Returns:
            flag: Variable that indicates to the user on return what action he has to do,
                or the state of the algorithm
        """
        lb_pass = lb.ctypes.data_as(POINTER(c_float)) if lb is not None else lb
        ub_pass = ub.ctypes.data_as(POINTER(c_float)) if ub is not None else ub
        cflag = c_int(flag)
        self._ctypes_lbfgs(
            n,
            x.ctypes.data_as(POINTER(c_float)),
            byref(c_float(fcost)),
            grad.ctypes.data_as(POINTER(c_float)),
            byref(self.udf),
            byref(cflag),
            lb_pass,
            ub_pass,
        )
        return cflag.value

    def TRN(
        self,
        n: int,
        x: ArrayLike,
        fcost: float,
        grad: ArrayLike,
        d: ArrayLike,
        hd: ArrayLike,
        flag: int,
        lb: Optional[ArrayLike] = None,
        ub: Optional[ArrayLike] = None,
    ) -> int:
        """Truncated-Newton algorithm.

        Args:
            n: The parameter space dimension (size of the x vector)
            x: Solution vector
            fcost: Cost function value associated with x
            grad: Gradient at x
            d: Vector used in the conjugate gradient process. At the beginning
                it is the negative residual, then it is updated after every iteration
                with the new search direction
            hd: Vector resulting from multiplication of the `d` vector and the
                Hessian operator
            flag: Variable that indicates to the user on return what action he has to do,
                or the state of the algorithm
            lb: Vector of lower bounds for x. By default None
            ub: Vector of upper bounds for x. By default None

        Returns:
            flag: Variable that indicates to the user on return what action he has to do,
                or the state of the algorithm
        """
        lb_pass = lb.ctypes.data_as(POINTER(c_float)) if lb is not None else lb
        ub_pass = ub.ctypes.data_as(POINTER(c_float)) if ub is not None else ub
        cflag = c_int(flag)
        self._ctypes_trn(
            n,
            x.ctypes.data_as(POINTER(c_float)),
            byref(c_float(fcost)),
            grad.ctypes.data_as(POINTER(c_float)),
            d.ctypes.data_as(POINTER(c_float)),
            hd.ctypes.data_as(POINTER(c_float)),
            byref(self.udf),
            byref(cflag),
            lb_pass,
            ub_pass,
        )
        return cflag.value

    def PTRN(
        self,
        n: int,
        x: ArrayLike,
        fcost: float,
        grad: ArrayLike,
        grad_preco: ArrayLike,
        residual: ArrayLike,
        residual_preco: ArrayLike,
        d: ArrayLike,
        hd: ArrayLike,
        flag: int,
        lb: Optional[ArrayLike] = None,
        ub: Optional[ArrayLike] = None,
    ) -> int:
        """Preconditioned Truncated-Newton algorithm.

        Args:
            n: The parameter space dimension (size of the x vector)
            x: Solution vector
            fcost: Cost function value associated with x
            grad: Gradient at x
            grad_preco: Preconditioned gradient at x
            residual: Residual vector used in the conjugate gradient process
            residual_preco: Preconditioned residual
            d: Vector used in the conjugate gradient process. At the beginning
                it is the negative preconditioned residual, then it is updated after
                every iteration with the new search direction
            hd: Vector resulting from multiplication of the `d` vector and the
                Hessian operator
            flag: Variable that indicates to the user on return what action he has to do,
                or the state of the algorithm
            lb: Vector of lower bounds for x. By default None
            ub: Vector of upper bounds for x. By default None

        Returns:
            flag: Variable that indicates to the user on return what action he has to do,
                or the state of the algorithm
        """
        lb_pass = lb.ctypes.data_as(POINTER(c_float)) if lb is not None else lb
        ub_pass = ub.ctypes.data_as(POINTER(c_float)) if ub is not None else ub
        cflag = c_int(flag)
        self._ctypes_ptrn(
            n,
            x.ctypes.data_as(POINTER(c_float)),
            byref(c_float(fcost)),
            grad.ctypes.data_as(POINTER(c_float)),
            grad_preco.ctypes.data_as(POINTER(c_float)),
            residual.ctypes.data_as(POINTER(c_float)),
            residual_preco.ctypes.data_as(POINTER(c_float)),
            d.ctypes.data_as(POINTER(c_float)),
            hd.ctypes.data_as(POINTER(c_float)),
            byref(self.udf),
            byref(cflag),
            lb_pass,
            ub_pass,
        )
        return cflag.value
