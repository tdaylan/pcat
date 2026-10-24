"""Point-spread-function utilities."""

import numpy as np
from copy import deepcopy


def _subpixel_cubic_design(x, y):
	x, y = np.broadcast_arrays(x, y)
	return np.moveaxis(np.array((np.ones_like(x), x, y, x*x, x*y, y*y,
	                            x*x*x, x*x*y, x*y*y, y*y*y)), 0, -1)


class CubicSubpixelPSF:
	"""Fit Cartesian PSF values as cubic functions of source subpixel shifts.

	The fixed radial callback returns (energy, position, quality) densities.
	Pixel coordinates lie on a regular lattice with the supplied pixel size.
	The ten-term fit is checked at interleaved phase samples. Cells exceeding
	0.1 percent of the profile peak use direct evaluation, as do positions
	outside the table and negative/nonfinite polynomial values. This sampled
	check is not a rigorous uniform error bound. No flux renormalization or
	value clipping is applied. Rebuild the object when PSF parameters change.
	"""

	def __init__(self, profile, x_coordinates, y_coordinates, pixel_size, relative_tolerance=1e-3):
		self._profile = deepcopy(profile)
		coordinates = []
		for values in (x_coordinates, y_coordinates):
			if np.ma.is_masked(values) or np.iscomplexobj(values):
				raise ValueError('Pixel coordinates must be unmasked real values')
			values = np.array(values, dtype=float, copy=True)
			if values.ndim != 1 or not values.size or not np.isfinite(values).all():
				raise ValueError('Pixel coordinates must be nonempty finite vectors')
			coordinates.append(values)
		if coordinates[0].shape != coordinates[1].shape:
			raise ValueError('Pixel coordinate vectors must have matching shapes')
		for value in (pixel_size, relative_tolerance):
			if (np.ma.is_masked(value) or np.iscomplexobj(value)
			        or np.ndim(value) != 0 or isinstance(value, (bool, np.bool_))
			        or not np.isfinite(value) or value <= 0):
				raise ValueError('Pixel size and relative tolerance must be positive finite scalars')
		self.pixel_size = float(pixel_size)
		self._coordinates = tuple(coordinates)
		self._origins = tuple(float(values.min()) for values in coordinates)
		self._pixel_indices = []
		for values, origin in zip(coordinates, self._origins):
			indices = (values - origin) / self.pixel_size
			if not np.allclose(indices, np.rint(indices), rtol=0, atol=1e-7):
				raise ValueError('Pixel coordinates must lie on the regular pixel lattice')
			self._pixel_indices.append(np.rint(indices).astype(int))
		self._spans = tuple(int(indices.max()) + 2 for indices in self._pixel_indices)
		offset_x, offset_y = np.meshgrid(
			np.arange(-self._spans[0], self._spans[0] + 1),
			np.arange(-self._spans[1], self._spans[1] + 1), indexing='ij',
		)
		self._table_shape = offset_x.shape
		phases = np.linspace(-0.5, 0.5, 5)
		phase_x, phase_y = np.meshgrid(phases, phases, indexing='ij')
		design = _subpixel_cubic_design(phase_x.ravel(), phase_y.ravel())

		def densities(phase_x, phase_y):
			distances = self.pixel_size * np.hypot(
				offset_x.ravel()[None, :] - phase_x[:, None],
				offset_y.ravel()[None, :] - phase_y[:, None],
			)
			values = self._profile(distances.ravel())
			if np.ma.is_masked(values) or np.iscomplexobj(values):
				raise ValueError('PSF callback must return unmasked real densities')
			values = np.asarray(values, dtype=float)
			if values.ndim != 3 or values.shape[1] != distances.size:
				raise ValueError('PSF callback must return energy-by-position-by-quality values')
			if not np.isfinite(values).all() or np.any(values < 0):
				raise ValueError('PSF callback must return finite nonnegative densities')
			return values.reshape(values.shape[0], len(phase_x), offset_x.size, values.shape[2]).transpose(1, 2, 0, 3)

		training = densities(phase_x.ravel(), phase_y.ravel())
		self._output_shape = training.shape[2:]
		coefficients = np.linalg.lstsq(design, training.reshape(25, -1), rcond=None)[0]
		self._coefficients = np.moveaxis(
			coefficients.reshape(10, *self._table_shape, *self._output_shape), 0, -1,
		).copy()
		validation = np.array([-0.375, -0.125, 0.125, 0.375])
		check_x, check_y = np.meshgrid(validation, validation, indexing='ij')
		truth = densities(check_x.ravel(), check_y.ravel())
		predicted = (_subpixel_cubic_design(check_x.ravel(), check_y.ravel()) @ coefficients).reshape(truth.shape)
		peak = np.max(training, axis=(0, 1))
		self._direct_cells = np.any(np.abs(predicted - truth) > relative_tolerance * peak, axis=(0, 2, 3)).reshape(self._table_shape)
		for array in (*coordinates, *self._pixel_indices, self._coefficients, self._direct_cells):
			array.setflags(write=False)

	def __deepcopy__(self, memo):
		return self

	def __call__(self, distances):
		return self._profile(distances)

	def evaluate(self, xpos, ypos, pixel_indices):
		if any(np.ma.is_masked(value) or np.iscomplexobj(value)
		       or np.ndim(value) != 0 or not np.isfinite(value) for value in (xpos, ypos)):
			raise ValueError('Source positions must be finite scalars')
		indices = np.asarray(pixel_indices)
		if indices.ndim != 1 or indices.dtype.kind not in 'iu' or np.any(indices < 0) or np.any(indices >= len(self._coordinates[0])):
			raise ValueError('Pixel indices must select valid grid positions')
		with np.errstate(over='ignore'):
			positions = (np.array([xpos, ypos]) - self._origins) / self.pixel_size
		if not np.isfinite(positions).all() or np.any(np.abs(positions) > np.iinfo(np.intp).max / 4):
			distances = np.hypot(self._coordinates[0][indices] - xpos, self._coordinates[1][indices] - ypos)
			return self._profile(distances)
		anchors = np.rint(positions).astype(int)
		phases = positions - anchors
		x_index = self._pixel_indices[0][indices] - anchors[0] + self._spans[0]
		y_index = self._pixel_indices[1][indices] - anchors[1] + self._spans[1]
		inside = (x_index >= 0) & (x_index < self._table_shape[0]) & (y_index >= 0) & (y_index < self._table_shape[1])
		result = np.zeros((self._output_shape[0], indices.size, self._output_shape[1]))
		if np.any(inside):
			selected = np.flatnonzero(inside)
			coefficient_rows = self._coefficients[x_index[inside], y_index[inside]]
			result[:, inside, :] = np.einsum('...k,k->...', coefficient_rows, _subpixel_cubic_design(*phases)).transpose(1, 0, 2)
			inside[selected[self._direct_cells[x_index[inside], y_index[inside]]]] = False
		inside &= np.isfinite(result).all(axis=(0, 2)) & (result >= 0).all(axis=(0, 2))
		if not np.all(inside):
			selected = indices[~inside]
			distances = np.hypot(self._coordinates[0][selected] - xpos, self._coordinates[1][selected] - ypos)
			result[:, ~inside, :] = self._profile(distances)
		return result


def psf_poly_fit(gdat, psfnusam, factusam):
	"""Fit a cubic subpixel model to each pixel of an oversampled 1D PSF.

	The returned columns contain coefficients in ascending polynomial order.
	Evaluating a column at offsets from zero through one reconstructs the
	oversampled PSF values across the corresponding detector pixel.
	"""
	psfnusam = np.asarray(psfnusam, dtype=np.float32)
	if psfnusam.ndim != 1:
		raise ValueError('psfnusam must be one-dimensional.')
	if not isinstance(factusam, (int, np.integer)) or factusam < 3:
		raise ValueError('factusam must be an integer of at least 3.')
	if psfnusam.size % factusam:
		raise ValueError('psfnusam.size must be divisible by factusam.')

	gdat.enerpsfnusam = np.arange(psfnusam.size + 1) - psfnusam.size / 2.
	psfnusampadd = np.pad(psfnusam, (0, 1))
	offs = np.linspace(0., 1., factusam + 1, dtype=np.float32)
	design = np.column_stack((np.ones_like(offs), offs, offs**2, offs**3))
	numbsidepsfn = psfnusam.size // factusam
	coefspix = np.empty((design.shape[1], numbsidepsfn), dtype=np.float32)

	for indxpixl in range(numbsidepsfn):
		strt = indxpixl * factusam
		stop = strt + factusam + 1
		coefspix[:, indxpixl] = np.linalg.lstsq(
			design, psfnusampadd[strt:stop], rcond=None,
		)[0]

	return coefspix