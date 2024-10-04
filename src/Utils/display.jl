"""Some convenience functions for displaying structures"""

# format floats
sfmt(f::AbstractFloat) = @sprintf("%.5g", f)

# fallback for sfmt
sfmt(x) = string(x)
