# ALRD News

## ALRD 0.0.0.9000

Initial development version.

### Package functions

- Added `simulate_alrd_data()` for toy longitudinal ITE data generation.
- Added `fit_lrd()` for longitudinal residual difference learning.
- Added `fit_alrd()` for augmented longitudinal residual difference learning.
- Added lasso and non-penalized ALRD fitting options.
- Added `predict_ite()` for ITE prediction from fitted LRD/ALRD objects.
- Added `estimate_ps()` for logistic propensity-score estimation.
- Added `working_covariance()` for inspecting subject-level working covariance
  matrices.

### Method support

- Supports treatment coded as `0/1` or `-1/1`.
- Supports independence, exchangeable, and AR(1) working covariance structures.
- Supports fixed `rho` or iterative residual-based `rho` estimation.
- Supports lasso selection with optional unpenalized refit.

### Documentation and development

- Added README quick-start examples.
- Added package usage guide.
- Added collaborator onboarding and GitHub workflow documentation.
- Added unit-test infrastructure.

