# This file was automatically generated by SWIG (http://www.swig.org).
# Version 4.0.2
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info < (2, 7, 0):
    raise RuntimeError("Python 2.7 or later required")

# Import the low-level C/C++ module
if __package__ or "." in __name__:
    from . import _distlink
else:
    import _distlink

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)


def _swig_setattr_nondynamic_instance_variable(set):
    def set_instance_attr(self, name, value):
        if name == "thisown":
            self.this.own(value)
        elif name == "this":
            set(self, name, value)
        elif hasattr(self, name) and isinstance(getattr(type(self), name), property):
            set(self, name, value)
        else:
            raise AttributeError("You cannot add instance attributes to %s" % self)
    return set_instance_attr


def _swig_setattr_nondynamic_class_variable(set):
    def set_class_attr(cls, name, value):
        if hasattr(cls, name) and not isinstance(getattr(cls, name), property):
            set(cls, name, value)
        else:
            raise AttributeError("You cannot add class attributes to %s" % cls)
    return set_class_attr


def _swig_add_metaclass(metaclass):
    """Class decorator for adding a metaclass to a SWIG wrapped class - a slimmed down version of six.add_metaclass"""
    def wrapper(cls):
        return metaclass(cls.__name__, cls.__bases__, cls.__dict__.copy())
    return wrapper


class _SwigNonDynamicMeta(type):
    """Meta class to enforce nondynamic attributes (no new attributes) for a class"""
    __setattr__ = _swig_setattr_nondynamic_class_variable(type.__setattr__)



def detect_suitable_options(max_root_error, min_root_error, max_anom_error):
    r"""detect_suitable_options(max_root_error, min_root_error, max_anom_error)"""
    return _distlink.detect_suitable_options(max_root_error, min_root_error, max_anom_error)
class COrbitData(object):
    r"""
    Proxy of C++ COrbitData< double > class.
    Proxy of C++ COrbitData< double > class.
    """

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def __init__(self, *args):
        r"""
        __init__(self) -> COrbitData
        __init__(self, a_, e_, i_, w_, Om_) -> COrbitData
        """
        _distlink.COrbitData_swiginit(self, _distlink.new_COrbitData(*args))

    def set_data(self, a_, e_, i_, w_, Om_):
        r"""set_data(self, a_, e_, i_, w_, Om_)"""
        return _distlink.COrbitData_set_data(self, a_, e_, i_, w_, Om_)

    def get_data(self, a_, e_, i_, w_, Om_):
        r"""get_data(self, a_, e_, i_, w_, Om_)"""
        return _distlink.COrbitData_get_data(self, a_, e_, i_, w_, Om_)

    def get_a(self):
        r"""get_a(self) -> double"""
        return _distlink.COrbitData_get_a(self)

    def get_e(self):
        r"""get_e(self) -> double"""
        return _distlink.COrbitData_get_e(self)

    def get_i(self):
        r"""get_i(self) -> double"""
        return _distlink.COrbitData_get_i(self)

    def get_w(self):
        r"""get_w(self) -> double"""
        return _distlink.COrbitData_get_w(self)

    def get_Om(self):
        r"""get_Om(self) -> double"""
        return _distlink.COrbitData_get_Om(self)

    def vectorP(self):
        r"""vectorP(self) -> double const *"""
        return _distlink.COrbitData_vectorP(self)

    def vectorQ(self):
        r"""vectorQ(self) -> double const *"""
        return _distlink.COrbitData_vectorQ(self)

    def get_vectors(self, P_, Q_):
        r"""get_vectors(self, P_, Q_)"""
        return _distlink.COrbitData_get_vectors(self, P_, Q_)
    __swig_destroy__ = _distlink.delete_COrbitData

# Register COrbitData in _distlink:
_distlink.COrbitData_swigregister(COrbitData)


def test_peri_apo(O1, O2, limit):
    r"""test_peri_apo(O1, O2, limit) -> bool"""
    return _distlink.test_peri_apo(O1, O2, limit)
class SMOIDResult(object):
    r"""
    Proxy of C++ SMOIDResult< double > class.
    Proxy of C++ SMOIDResult< double > class.
    """

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    good = property(_distlink.SMOIDResult_good_get, _distlink.SMOIDResult_good_set, doc=r"""good""")
    distance = property(_distlink.SMOIDResult_distance_get, _distlink.SMOIDResult_distance_set, doc=r"""distance""")
    distance_error = property(_distlink.SMOIDResult_distance_error_get, _distlink.SMOIDResult_distance_error_set, doc=r"""distance_error""")
    u1 = property(_distlink.SMOIDResult_u1_get, _distlink.SMOIDResult_u1_set, doc=r"""u1""")
    u1_error = property(_distlink.SMOIDResult_u1_error_get, _distlink.SMOIDResult_u1_error_set, doc=r"""u1_error""")
    u2 = property(_distlink.SMOIDResult_u2_get, _distlink.SMOIDResult_u2_set, doc=r"""u2""")
    u2_error = property(_distlink.SMOIDResult_u2_error_get, _distlink.SMOIDResult_u2_error_set, doc=r"""u2_error""")
    root_count = property(_distlink.SMOIDResult_root_count_get, _distlink.SMOIDResult_root_count_set, doc=r"""root_count""")
    min_delta = property(_distlink.SMOIDResult_min_delta_get, _distlink.SMOIDResult_min_delta_set, doc=r"""min_delta""")
    iter_count = property(_distlink.SMOIDResult_iter_count_get, _distlink.SMOIDResult_iter_count_set, doc=r"""iter_count""")
    iter_count_2D = property(_distlink.SMOIDResult_iter_count_2D_get, _distlink.SMOIDResult_iter_count_2D_set, doc=r"""iter_count_2D""")
    time = property(_distlink.SMOIDResult_time_get, _distlink.SMOIDResult_time_set, doc=r"""time""")

    def __init__(self):
        r"""__init__(self) -> SMOIDResult"""
        _distlink.SMOIDResult_swiginit(self, _distlink.new_SMOIDResult())
    __swig_destroy__ = _distlink.delete_SMOIDResult

# Register SMOIDResult in _distlink:
_distlink.SMOIDResult_swigregister(SMOIDResult)

class SLCResult(object):
    r"""
    Proxy of C++ SLCResult< double > class.
    Proxy of C++ SLCResult< double > class.
    """

    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    I = property(_distlink.SLCResult_I_get, _distlink.SLCResult_I_set, doc=r"""I""")
    l = property(_distlink.SLCResult_l_get, _distlink.SLCResult_l_set, doc=r"""l""")
    lmod = property(_distlink.SLCResult_lmod_get, _distlink.SLCResult_lmod_set, doc=r"""lmod""")
    l2 = property(_distlink.SLCResult_l2_get, _distlink.SLCResult_l2_set, doc=r"""l2""")

    def __init__(self):
        r"""__init__(self) -> SLCResult"""
        _distlink.SLCResult_swiginit(self, _distlink.new_SLCResult())
    __swig_destroy__ = _distlink.delete_SLCResult

# Register SLCResult in _distlink:
_distlink.SLCResult_swigregister(SLCResult)


def LC(O1, O2, min_mut_incl):
    r"""LC(O1, O2, min_mut_incl) -> SLCResult"""
    return _distlink.LC(O1, O2, min_mut_incl)

def MOID_fast(*args):
    r"""MOID_fast(O1, O2, max_root_error, min_root_error, nu=static_cast< double >(1)) -> SMOIDResult"""
    return _distlink.MOID_fast(*args)

def MOID_direct_search(O1, O2, densities, max_dist_error, max_anom_error):
    r"""MOID_direct_search(O1, O2, densities, max_dist_error, max_anom_error) -> SMOIDResult"""
    return _distlink.MOID_direct_search(O1, O2, densities, max_dist_error, max_anom_error)

def restrict_search_range(O1, O2, u1, u1_):
    r"""restrict_search_range(O1, O2, u1, u1_) -> int"""
    return _distlink.restrict_search_range(O1, O2, u1, u1_)


