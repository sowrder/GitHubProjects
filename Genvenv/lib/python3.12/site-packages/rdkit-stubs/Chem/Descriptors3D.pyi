"""
 Descriptors derived from a molecule's 3D structure

"""
from __future__ import annotations
from rdkit.Chem.Descriptors import _isCallable
from rdkit.Chem import rdMolDescriptors
__all__ = ['CalcMolDescriptors3D', 'descList', 'rdMolDescriptors']
def CalcMolDescriptors3D(mol, confId = None):
    """
    
        Compute all 3D descriptors of a molecule
        
        Arguments:
        - mol: the molecule to work with
        - confId: conformer ID to work with. If not specified the default (-1) is used
        
        Return:
        
        dict
            A dictionary with decriptor names as keys and the descriptor values as values
    
        raises a ValueError 
            If the molecule does not have conformers
        
    """
def _setupDescriptors(namespace):
    ...
descList: list  # value = [('PMI1', <function <lambda> at 0x104b32340>), ('PMI2', <function <lambda> at 0x107a96520>), ('PMI3', <function <lambda> at 0x107a965c0>), ('NPR1', <function <lambda> at 0x107a96660>), ('NPR2', <function <lambda> at 0x107a96700>), ('RadiusOfGyration', <function <lambda> at 0x107a967a0>), ('InertialShapeFactor', <function <lambda> at 0x107a96840>), ('Eccentricity', <function <lambda> at 0x107a968e0>), ('Asphericity', <function <lambda> at 0x107a96980>), ('SpherocityIndex', <function <lambda> at 0x107a96a20>), ('PBF', <function <lambda> at 0x107a96ac0>)]
