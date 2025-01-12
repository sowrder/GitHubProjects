from __future__ import annotations
import typing
__all__ = ['AlignMol', 'FloatVector', 'PrepareConformer', 'ShapeInput']
class FloatVector(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
            C++ signature :
                bool __contains__(std::__1::vector<float, std::__1::allocator<float>> {lvalue},_object*)
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
            C++ signature :
                void __delitem__(std::__1::vector<float, std::__1::allocator<float>> {lvalue},_object*)
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
            C++ signature :
                boost::python::api::object __getitem__(boost::python::back_reference<std::__1::vector<float, std::__1::allocator<float>>&>,_object*)
        """
    def __init__(self) -> None:
        """
            C++ signature :
                void __init__(_object*)
        """
    def __iter__(self) -> typing.Any:
        """
            C++ signature :
                boost::python::objects::iterator_range<boost::python::return_value_policy<boost::python::return_by_value, boost::python::default_call_policies>, std::__1::__wrap_iter<float*>> __iter__(boost::python::back_reference<std::__1::vector<float, std::__1::allocator<float>>&>)
        """
    def __len__(self) -> int:
        """
            C++ signature :
                unsigned long __len__(std::__1::vector<float, std::__1::allocator<float>> {lvalue})
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
            C++ signature :
                void __setitem__(std::__1::vector<float, std::__1::allocator<float>> {lvalue},_object*,_object*)
        """
    def append(self, item: typing.Any) -> None:
        """
            C++ signature :
                void append(std::__1::vector<float, std::__1::allocator<float>> {lvalue},boost::python::api::object)
        """
    def extend(self, other: typing.Any) -> None:
        """
            C++ signature :
                void extend(std::__1::vector<float, std::__1::allocator<float>> {lvalue},boost::python::api::object)
        """
class ShapeInput(Boost.Python.instance):
    @staticmethod
    def __init__(*args, **kwargs):
        """
        Raises an exception
        This class cannot be instantiated from Python
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @property
    def alpha_vector(*args, **kwargs):
        ...
    @alpha_vector.setter
    def alpha_vector(*args, **kwargs):
        ...
    @property
    def atom_type_vector(*args, **kwargs):
        ...
    @atom_type_vector.setter
    def atom_type_vector(*args, **kwargs):
        ...
    @property
    def coord(*args, **kwargs):
        ...
    @coord.setter
    def coord(*args, **kwargs):
        ...
    @property
    def shift(*args, **kwargs):
        ...
    @shift.setter
    def shift(*args, **kwargs):
        ...
    @property
    def sof(*args, **kwargs):
        ...
    @sof.setter
    def sof(*args, **kwargs):
        ...
    @property
    def sov(*args, **kwargs):
        ...
    @sov.setter
    def sov(*args, **kwargs):
        ...
    @property
    def volumeAtomIndexVector(*args, **kwargs):
        ...
    @volumeAtomIndexVector.setter
    def volumeAtomIndexVector(*args, **kwargs):
        ...
@typing.overload
def AlignMol(ref: Mol, probe: Mol, refConfId: int = -1, probeConfId: int = -1, useColors: bool = True, opt_param: float = 0.5, max_preiters: int = 3, max_postiters: int = 16) -> tuple:
    """
        aligns probe to ref, probe is modified
    
        C++ signature :
            boost::python::tuple AlignMol(RDKit::ROMol,RDKit::ROMol {lvalue} [,int=-1 [,int=-1 [,bool=True [,double=0.5 [,unsigned int=3 [,unsigned int=16]]]]]])
    """
@typing.overload
def AlignMol(refShape: ShapeInput, probe: Mol, probeConfId: int = -1, useColors: bool = True, opt_param: float = 0.5, max_preiters: int = 3, max_postiters: int = 16) -> tuple:
    """
        aligns probe to reference shape, probe is modified
    
        C++ signature :
            boost::python::tuple AlignMol(ShapeInput,RDKit::ROMol {lvalue} [,int=-1 [,bool=True [,double=0.5 [,unsigned int=3 [,unsigned int=16]]]]])
    """
def PrepareConformer(mol: Mol, confId: int = -1, useColors: bool = True) -> ShapeInput:
    """
        returns a shape object for a molecule
    
        C++ signature :
            ShapeInput* PrepareConformer(RDKit::ROMol [,int=-1 [,bool=True]])
    """
