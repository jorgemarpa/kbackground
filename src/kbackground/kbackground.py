from scipy import sparse
import numpy as np
from patsy import dmatrix
from dataclasses import dataclass
from astropy.io import fits
from astropy.stats import sigma_clip


@dataclass
class Estimator:
    """Background Estimator for Kepler/K2

    Parameters
    ----------


    """
    row:np.ndarray
    column:np.ndarray
    flux:np.ndarray
    cadenceno:np.ndarray

    def __post_init__(self):
        self.xknots, self.yknots = xknots, yknots = np.linspace(20, 1108, 42)[1:-1], np.linspace(27, 1040, 42)[1:-1]
        self.mask = ~sigma_clip(np.median(self.flux, axis=0)).mask

        med_flux = np.median(self.flux, axis=0)[None, :]
        self.flux_offset = np.median(self.flux - med_flux, axis=1)

        self.A = self._make_A(self.row, self.column)
        prior_mu = np.zeros(self.A.shape[1])
        prior_mu[0] = 1
        prior_mu = (self.flux_offset[:, None] * prior_mu)
        prior_sigma = np.ones(self.A.shape[1]) * 40

        self.sigma_w_inv = self.A[self.mask].T.dot(self.A[self.mask]) + np.diag(1 / prior_sigma ** 2)
        Bs = self.A[self.mask].T.dot((self.flux - med_flux)[:, self.mask].T) + (prior_mu/prior_sigma**2).T
        self.ws = np.linalg.solve(self.sigma_w_inv, Bs).T
        self._model_row = self.row
        self._model_column = self.column
        self._model_A = self.A

    @staticmethod
    def from_mission_bkg(fname):
        hdu = fits.open(fname)
        self = Estimator(hdu[2].data['RAWX'], hdu[2].data['RAWY'], hdu[1].data['FLUX'], hdu[1].data['CADENCENO'])
        hdr = hdu[0].header
        self.channel = hdr['CHANNEL']
        self.mission = hdr['MISSION']
        if 'QUARTER' in hdr:
            self.quarter = hdr['QUARTER']
        if 'CAMPAIGN' in hdr:
            self.campaign = hdr['CAMPAIGN']
        return self

    def model(self, cadenceno, row=None, column=None):
        if row is not None:
            if (self._model_row is None) | np.atleast_1d(((self._model_row != row) | (self._model_column != column))).any():
                self._model_row = row
                self._model_column = column
                self._model_A = self._make_A(row, column)
            return self._model_A.dot(self.ws[np.in1d(self.cadenceno, cadenceno)][0])
        else:
            return self.A.dot(self.ws[np.in1d(self.cadenceno, cadenceno)][0])

    def __repr__(self):
        if hasattr(self, 'quarter'):
            return f'KBackground.Estimator Channel:{self.channel} Quarter:{self.quarter}'
        if hasattr(self, 'campaign'):
            return f'KBackground.Estimator Channel:{self.channel} Campaign:{self.campaign}'
        return f'KBackground.Estimator'

    @property
    def shape(self):
        return self.flux.shape

    def _make_A(self, x, y):
        x_spline = sparse.csr_matrix(
            np.asarray(
                dmatrix(
                    "bs(x, knots=knots, degree=3, include_intercept=True)",
                    {"x": np.hstack([0, x, 1400]), "knots": self.xknots},
                )
            )
        )[1:-1]

        # y_spline = sparse.csr_matrix(
        #     np.asarray(
        #         dmatrix(
        #             "bs(x, knots=knots, degree=3, include_intercept=True)",
        #             {"x": np.hstack([0, list(y), 1400]), "knots": self.yknots},
        #         )
        #     )
        # )[1:-1]
        # X = sparse.hstack(
        #     [x_spline.multiply(y_spline[:, idx]) for idx in range(y_spline.shape[1])],
        #     format="csr",
        # )
        X = sparse.hstack(
            [x_spline,
             #y_spline,
            # #x_spline.multiply(y[:, None] - y.mean()),
            # y_spline.multiply(x[:, None] - x.mean()),
            # y_spline.multiply((x[:, None] - x.mean())**2),
            ],
            format="csr",
        )

        return X
