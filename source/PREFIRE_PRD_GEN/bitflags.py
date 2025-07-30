"""
Module containing methods to work with "bitflags", where the individual binary
 bits in an integer datatype contain boolean info.

This program requires python version 3.6 or later, and is importable as a 
python module.
"""

#--------------------------------------------------------------------------
def bit_meaning(int_v, offset):
   """Returns True if the bit at 'offset' is one, and False otherwise."""
   mask = 1 << offset
   if (int_v & mask) != 0:  # nonzero result (2**offset) if the bit is one
      return True
   else:
      return False


#--------------------------------------------------------------------------
def set_bit(int_v, offset):
   """Returns an integer with the bit at 'offset' set to one."""
   mask = 1 << offset
   return(int_v | mask)


#--------------------------------------------------------------------------
def clear_bit(int_v, offset):
    """Returns an integer with the bit at 'offset' cleared (set to zero)."""
    mask = ~(1 << offset)
    return (int_v & mask)


#--------------------------------------------------------------------------
def toggle_bit(int_v, offset):
    """Returns an integer with the bit at 'offset' inverted, 0 -> 1, 1 -> 0."""
    mask = 1 << offset
    return (int_v ^ mask)
