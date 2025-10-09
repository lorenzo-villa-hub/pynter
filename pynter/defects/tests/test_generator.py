
from pymatgen.io.vasp.inputs import Poscar

from pynter.testing.core import PynterTest
from pynter.defects.generator import create_vacancies, create_substitutions


class TestDefectGenerator(PynterTest):

    @classmethod
    def setUpClass(cls):
        self = cls()
        structure = Poscar.from_file(self.get_testfile_path('BaTiO3_mp-5986.poscar')).structure
        cls.structure = structure

    
    def test_create_vacancies(self):
        defects = create_vacancies(
                                structure=self.structure,
                                elements=None,
                                supercell_size=3)
        
        self.assertEqual( len(defects), 4)
        
        def count_species(defects,element):
            count = 0
            for df in defects:
                if df.specie == element:
                    count += 1
            return count
        
        self.assertEqual(count_species(defects,'Ba'), 1)
        self.assertEqual(count_species(defects,'Ti'), 1)
        self.assertEqual(count_species(defects,'O'), 2)



    def test_create_substitutions(self):
        defects = create_substitutions(
                                structure=self.structure,
                                elements_to_replace={'Ba':'Sr','Ti':'Zr','O':'N'},
                                supercell_size=3)
        
        self.assertEqual( len(defects), 4)
        
        def count_species(defects,element):
            count = 0
            for df in defects:
                if df.specie == element:
                    count += 1
            return count
        
        self.assertEqual(count_species(defects,'Sr'), 1)
        self.assertEqual(count_species(defects,'Zr'), 1)
        self.assertEqual(count_species(defects,'N'), 2)