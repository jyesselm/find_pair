# Validation Package Usage

The validation package provides clean, focused modules for geometric validation of base pairs.

## Quick Start

```python
from prototypes.pair_identification.validation import (
    GeometricValidator,
    ValidationThresholds,
    ValidationResult,
)

# Use default thresholds
validator = GeometricValidator()
result = validator.validate(frame1, frame2, n1_pos, n9_pos)

if result.is_valid:
    print(f"Valid pair with quality score: {result.quality_score:.2f}")
else:
    print(f"Invalid: {result.rejection_reason}")
    print(f"Failed checks: {result.get_failed_checks()}")
```

## Using Custom Thresholds

```python
# Use strict validation
strict_validator = GeometricValidator(ValidationThresholds.strict())

# Use relaxed validation
relaxed_validator = GeometricValidator(ValidationThresholds.relaxed())

# Custom thresholds
custom = ValidationThresholds(
    max_dorg=10.0,
    max_d_v=2.0,
    max_plane_angle=50.0,
    min_dNN=5.0,
)
custom_validator = GeometricValidator(custom)
```

## Threshold Presets

| Preset | max_dorg | max_d_v | max_plane_angle |
|--------|----------|---------|-----------------|
| strict | 12.0 | 2.0 | 45.0 |
| default | 15.0 | 2.5 | 65.0 |
| relaxed | 18.0 | 3.0 | 75.0 |

## ValidationResult Fields

```python
result = validator.validate(frame1, frame2, n1_pos, n9_pos)

# Geometric metrics
print(f"Origin distance: {result.dorg:.2f} Å")
print(f"Vertical distance: {result.d_v:.2f} Å")
print(f"Plane angle: {result.plane_angle:.1f}°")
print(f"N1/N9 distance: {result.dNN:.2f} Å")

# Direction cosines
print(f"X-axis alignment: {result.dir_x:.3f}")
print(f"Y-axis alignment: {result.dir_y:.3f}")
print(f"Z-axis alignment: {result.dir_z:.3f}")

# Quality and validation
print(f"Quality score: {result.quality_score:.3f}")
print(f"Valid: {result.is_valid}")

# Individual checks
print(f"Distance check: {result.distance_check}")
print(f"d_v check: {result.d_v_check}")
print(f"Plane angle check: {result.plane_angle_check}")
print(f"dNN check: {result.dNN_check}")
```

## Debugging Failed Validations

```python
result = validator.validate(frame1, frame2, n1_pos, n9_pos)

if not result.is_valid:
    # Get list of failed checks
    failed = result.get_failed_checks()
    print(f"Failed: {', '.join(failed)}")
    
    # Or use the pre-formatted reason
    print(result.rejection_reason)
    
    # Check specific failures
    if "dorg" in failed:
        print(f"Origin distance {result.dorg:.2f} exceeds threshold")
    if "d_v" in failed:
        print(f"Vertical distance {result.d_v:.2f} exceeds threshold")
```

## Module Structure

- `validation/result.py` - ValidationResult dataclass
- `validation/thresholds.py` - ValidationThresholds with presets
- `validation/validator.py` - GeometricValidator class
- `validation/__init__.py` - Public exports

All modules follow strict code quality guidelines:
- Max 3 levels of indentation
- All functions <= 30 lines
- Comprehensive docstrings with type hints
