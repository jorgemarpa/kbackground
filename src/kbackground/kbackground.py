from dataclasses import dataclass

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip
from patsy import dmatrix
from scipy import sparse


@dataclass
class Estimator:
    """Background Estimator for Kepler/K2

    Parameters
    ----------

    row: np.ndarray
        1D array of row positions for pixels to calculate the background model of with shape npixels
    column: np.ndarray
        1D array of column positions for pixels to calculate the background model of with shape npixels
    flux : np.ndarray
        2D array of fluxes with shape ntimes x npixels
    """

    row: np.ndarray
    column: np.ndarray
    flux: np.ndarray

    def __post_init__(self):
        self.xknots, self.yknots = (
            np.linspace(20, 1108, 42)[1:-1],
            np.linspace(27, 1040, 8)[1:-1],
        )
        med_flux = np.median(self.flux, axis=0)[None, :]
        f = self.flux - med_flux
        # Mask out pixels that are particularly bright.
        self.mask = (med_flux[0] - np.percentile(med_flux, 20)) < 30
        if not self.mask.any():
            raise ValueError("All the input pixels are brighter than 30 counts.")
        self.mask &= ~sigma_clip(med_flux[0]).mask
        self.mask &= ~sigma_clip(np.std(f, axis=0)).mask

        self.flux_offset = np.median(f, axis=1)

        self.A = self._make_A(self.row, self.column)
        prior_mu = np.zeros(self.A.shape[1])
        prior_mu[0] = 1
        prior_mu = self.flux_offset[:, None] * prior_mu
        # Hard coding a prior with 100 count width.
        prior_sigma = np.ones(self.A.shape[1]) * 100

        self.sigma_w_inv = self.A[self.mask].T.dot(self.A[self.mask]) + np.diag(
            1 / prior_sigma ** 2
        )
        Bs = (
            self.A[self.mask].T.dot((f)[:, self.mask].T)
            + (prior_mu / prior_sigma ** 2).T
        )
        self.ws = np.linalg.solve(self.sigma_w_inv, Bs)
        self._model_row = self.row
        self._model_column = self.column
        self._model_A = self.A

    @staticmethod
    def from_mission_bkg(fname):
        hdu = fits.open(fname)
        self = Estimator(
            hdu[2].data["RAWX"],
            hdu[2].data["RAWY"],
            hdu[1].data["FLUX"],
        )
        return self

    def model(self, index=None, row=None, column=None):
        """returns the background model

        Parameters
        ----------
        index : int
            Index to provide model. If None, will model at all provided cadences.
        row: np.ndarray
            The row to provide the model at. Must be 1D. If none, will use the training dataset.
        column: np.ndarray
            The column to provide the model at. Must be 1D. If none, will use the training dataset.

        Returns
        -------
        model: np.ndarray
            The model flux. 2D array with shape nindex x npixels.
        """
        if index is None:
            index = np.arange(self.shape[0])
        index = np.atleast_1d(index)
        if row is not None:
            if (self._model_row is None) | np.atleast_1d(
                ((self._model_row != row) | (self._model_column != column))
            ).any():
                self._model_row = row
                self._model_column = column
                self._model_A = self._make_A(row, column)
            return np.atleast_2d(self._model_A.dot(self.ws[:, index])).T
        else:
            return np.atleast_2d(self.A.dot(self.ws[:, index])).T

    def __repr__(self):
        return "KBackground.Estimator"

    @property
    def shape(self):
        return self.flux.shape

    def _make_A(self, x, y):
        """Makes a reasonable design matrix for the rolling band."""
        x_spline = sparse.csr_matrix(
            np.asarray(
                dmatrix(
                    "bs(x, knots=knots, degree=3, include_intercept=True)",
                    {"x": np.hstack([0, x, 1400]), "knots": self.xknots},
                )
            )
        )[1:-1]

        y_spline = sparse.csr_matrix(
            np.asarray(
                dmatrix(
                    "bs(x, knots=knots, degree=3, include_intercept=True)",
                    {"x": np.hstack([0, list(y), 1400]), "knots": self.yknots},
                )
            )
        )[1:-1]
        X = sparse.hstack(
            [x_spline.multiply(y_spline[:, idx]) for idx in range(y_spline.shape[1])],
            format="csr",
        )

        return X
