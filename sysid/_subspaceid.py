# -*- coding: utf-8 -*-
"""
Implementation of stochastic subspace identification methods
"""
import abc
import functools
import numpy as np
from .utils import find_psd_matrix, get_frequency_vector


__all__ = ["CovarianceDrivenStochasticSID",
           "AbstractReferenceBasedStochasticSID"]


def create_block_hankel_matrix(data, block_rows, ncols=None, ix_ref=None):
    """Block hankel matrix from data array

    Arguments
    ---------
    data : 2darray
        Data array where each row contains data from one
        sensor and each column corresponds to a specific
        time.
    block_rows : int
        Number of block rows
    ncols : int, optional
        Number of columns in block hankel matrix. If None,
        all data in the data matrix is used.
    ix_ref : list, optional
        Indices to the reference outputs in y. If `None`, all outputs
        are considered to be references.

    Returns
    -------
    2darray
        Block hankel matrix from data array
    """
    l, s = data.shape
    ix_ref = ix_ref or [*range(l)]
    r = len(ix_ref)
    i = block_rows
    j = ncols or s - 2*i + 1
    y = data
    yref = y[ix_ref]
    H = np.zeros(((r+l)*i, j))
    for m in range(2*i):
        if m < i:
            H[m*r:(m+1)*r, :] = yref[:, m:m+j]
        else:
            H[r*i+(m-i)*l:r*i+(m+1-i)*l, :] = y[:, m:m+j]
    return 1./np.sqrt(j)*H


class AbstractReferenceBasedStochasticSID(abc.ABC):
    @abc.abstractmethod
    def __init__(self, y, fs, ix_references=None):
        """Subspace identificator

        Arguments
        ---------
        y : 2darray
            Output data matrix (l x s) from `l` outputs with `s` samples.
        fs : float
            Sampling rate
        ix_references : list, optional
            Indices to the reference outputs in y. If `None`, all outputs
            are considered to be references.
        """
        self.y = y
        self.fs = fs
        self.ix_references = ix_references or [*range(self.l)]

    @abc.abstractmethod
    def perform(self, *args, **kwargs):
        pass

    def psdy(self, return_trace=False, **kw):
        """Compute power spectral density matrix of outputs

        Compute the power spectral density matrix of the outputs
        with Welch's method.

        Arguments
        ---------
        return_trace : bool, optional
            Return the entire psd matrix or only the trace of the psd matrix.
        kw : dict
            See keywords to scipy.signal.csd and scipy.signal.welch.

        Returns
        -------
        f : 1darray
            Frequency vector of the psd matrix
        Pyy : 3darray or 1darray
            Output PSD matrix where the first dimension refers to the
            frequency of the psd estimator, see get_frequency_vector,
            and the second and third dimensions refers to the degree
            of freedom of the input and output as given in y. If
            return_trace, the trace of the psd matrix is returned
            instead of the entire psd matrix.
        """
        psd = find_psd_matrix(self.y, self.y, fs=self.fs, **kw)
        f = get_frequency_vector(self.fs, psd.shape[2])
        if return_trace:
            out = np.trace(psd)
        else:
            out = psd
        return f, out

    @property
    def yref(self):
        return self.y[self.ix_references]

    @property
    def l(self):
        return self.y.shape[0]

    @property
    def r(self):
        return len(self.ix_references)

    @property
    def s(self):
        return self.y.shape[1]

    def j(self, i):
        """Number of columns in block hankel matrices

        Returns the number of columns in the block hankel
        matrices given that all the data is used in the
        identification problem

        Arguments
        ---------
        i : int
            Number of block rows

        Returns
        -------
        int
            Number of columns in the block hankel
            matrices.
        """
        return self.s-2*i+1

    @functools.lru_cache(maxsize=20, typed=False)
    def _Y(self, i):
        """Output block hankel matrix

        Arguments
        ---------
        i : int
            Number of block rows

        Returns
        -------
        2darray
            Output block hankel matrix
        """
        return create_block_hankel_matrix(
            self.y, i, ix_ref=self.ix_references)


class CovarianceDrivenStochasticSID(AbstractReferenceBasedStochasticSID):
    """Stochastic subspace identifier (SSI)

    Given measurements of output y, identify
    the system matrices A,C, and
    the covariance matrices Q, R, S of the process noise w and measurement noise v for the system.

        x_{k+1} = Ax_k + w_k
        y_k = Cx_k + v_k

    Implementation is based on [Overschee1996] and [Peeters1999].

    References
    ----------
    [Overschee1996] Van Overschee, P., De Moor, B., 1996.
        Subspace Identification for Linear Systems.
        Springer US, Boston, MA. doi: 10.1007/978-1-4613-0465-4

    [Peeters1999] Peeters, B., De Roeck, G., 1999.
        Reference based stochastic subspace identification
        for output only modal analysis.
        Mechanical Systems and Signal Processing 13, 855–878.
        doi: 10.1006/mssp.1999.1249
    """
    def __init__(self, y, fs, ix_references=None):
        """Reference-based covariance-driven stochastic subspace identifier.

        Define a reference-based covariance driven stochastic subspace identifier (SSI-cov/ref).
        See [Overschee1996] and [Peeters1999] for more information.

        Arguments
        ---------
        y : 2darray
            Output data matrix (l x s) from `l` outputs with `s` samples.
        fs : float
            Sampling rate (Hz)
        ix_references : list, optional
            Indices to the reference outputs in y. If `None`, all outputs
            are considered to be references.
        """
        self.y = y
        self.fs = fs
        self.ix_references = ix_references or [*range(self.l)] # l is the # of the nodes


    @functools.lru_cache(maxsize=20, typed=False)
    def _R(self, lag):
        """Correlation matrix of data array

        Correlation matrix between data array with `lag` last samples
        removed and data array with first `lag` samples removed.

        Arguments
        ---------
        lag : int
            Time lag / shift

        Returns
        -------
        2darray
            Correlation matrix
        """
        s = self.s
        i = np.abs(lag)
        return self.y[:, :s-i].dot(self.yref[:, i:].T) / (s-i)

    @functools.lru_cache(maxsize=20, typed=False)
    def _T(self, i):
        """Block toeplitz matrix from output correlations

        Arguments
        ---------
        i : int
            Number of block rows

        Returns
        -------
        2darray
            Block toeplitz matrix from output correlations
        """
        Y = self._Y(i)
        Yp = Y[:self.r*i]
        Yf = Y[self.r*i:]
        return Yf @ Yp.T

    @functools.lru_cache(maxsize=20, typed=False)
    def _svd_block_toeplitz(self, i):
        """Perform and return SVD of the block toeplitz matrix

        This method is cached, repeated calls with the same argument
        does not recompute the block toeplitz matrix.

        Arguments
        ---------
        block_rows : int
            Number of block rows

        Returns
        -------
        U, S, V : 2darray
            U and V are the unitary matrices and S is the
            singular values.
        """
        U, s, VH = np.linalg.svd(self._T(i))
        return U, s, VH

    def perform(self, order, block_rows):
        """Perform system identification

        Arguments
        ---------
        order : int
            Order of the identified model
        block_rows : int
            Number of block rows

        Returns
        -------
        A, C, G, R0 : 2darrays
            System, output, next state-output covariance
            and zero lag correlation matrix.

        Raises
        ------
        ValueError
            The ratio between the order and number of block rows must
            be less or equal than the number of references to ensure a
            consistent equation system, i.e. the system has valid
            dimensions.
        """
        i, l, n, r = block_rows, self.l, order, self.r
        if n/i > r:
            raise ValueError(
                "Following condition violated: order / block_rows <= r")
        U, s, VH = self._svd_block_toeplitz(i)
        U1 = U[:, :n]
        V1H = VH[:n]
        sqrt_S1 = np.diag(np.sqrt(s[:n]))
        inv_sqrt_S1 = np.diag(1/np.sqrt(s[:n]))

        Oi = U1 @ sqrt_S1
        C = Oi[:l]

        Cref = sqrt_S1 @ V1H
        G = Cref[:, -r:]

        T1 = self._T(i)
        T2 = np.zeros_like(T1)
        T2[:-l, :] = T1[l:, :]
        T2[-l:, r:] = T1[-l:, :-r]
        T2[-l:, :r] = self._R(2*i)
        A = inv_sqrt_S1 @ U1.T @ T2 @ V1H.T @ inv_sqrt_S1

        R0 = self.y @ self.y.T / self.s
        return A, C, G, R0
