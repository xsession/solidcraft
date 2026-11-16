#!/usr/bin/env python3
"""Run real cases and independently validate their XML arrays and field associations.

XML checks are always used. --vtk-reader additionally exercises the installed
VTK XML reader; neither mode claims a ParaView application/filter qualification.
"""
import argparse
import json
import math
from pathlib import Path
import subprocess
import tempfile
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]

def validate_vtu(path, vtk_reader=False):
    tree = ET.parse(path)
    assert tree.getroot().attrib['type'] == 'UnstructuredGrid'
    piece = tree.find('UnstructuredGrid/Piece')
    points = int(piece.attrib['NumberOfPoints'])
    cells = int(piece.attrib['NumberOfCells'])
    xyz = list(map(float, piece.find('Points/DataArray').text.split()))
    assert len(xyz) == points*3 and all(map(math.isfinite, xyz))
    for kind, count in [('PointData', points), ('CellData', cells)]:
        for array in piece.findall(kind+'/DataArray'):
            values = list(map(float, (array.text or '').split()))
            assert len(values) == count*int(array.attrib.get('NumberOfComponents', 1))
            assert all(map(math.isfinite, values))
    arrays = {a.attrib['Name']: list(map(int, (a.text or '').split())) for a in piece.findall('Cells/DataArray')}
    assert len(arrays['offsets']) == len(arrays['types']) == cells
    assert all(0 <= n < points for n in arrays['connectivity'])
    arity = {1:1,3:2,5:3,9:4,10:4,12:8,13:6,14:5}
    previous = 0
    for offset, kind in zip(arrays['offsets'], arrays['types']):
        assert offset-previous == arity[kind]
        previous = offset
    assert previous == len(arrays['connectivity'])
    if vtk_reader:
        from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader
        reader = vtkXMLUnstructuredGridReader()
        errors = []
        reader.AddObserver('ErrorEvent', lambda *_: errors.append('VTK reader error'))
        reader.SetFileName(str(path))
        reader.Update()
        assert not errors, (path, errors)
        grid = reader.GetOutput()
        assert grid.GetNumberOfPoints() == points and grid.GetNumberOfCells() == cells
        for i in range(points):
            assert all(math.isclose(x, y, rel_tol=1e-13, abs_tol=1e-14)
                       for x, y in zip(grid.GetPoint(i), xyz[3*i:3*i+3]))
        for i, kind in enumerate(arrays['types']):
            assert grid.GetCellType(i) == kind
        for tag, data, count in [('PointData', grid.GetPointData(), points),
                                  ('CellData', grid.GetCellData(), cells)]:
            expected = piece.findall(tag+'/DataArray')
            assert data.GetNumberOfArrays() == len(expected)
            for xml_array in expected:
                array = data.GetArray(xml_array.attrib['Name'])
                components = int(xml_array.attrib.get('NumberOfComponents', 1))
                assert array is not None and array.GetNumberOfTuples() == count
                assert array.GetNumberOfComponents() == components
                values = list(map(float, (xml_array.text or '').split()))
                for index, value in enumerate(values):
                    assert math.isclose(array.GetComponent(index//components, index%components),
                                        value, rel_tol=1e-13, abs_tol=1e-14)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--keep', type=Path)
    p.add_argument('--vtk-reader', action='store_true')
    a = p.parse_args()
    binary = a.binary.resolve()
    checked = 0
    with tempfile.TemporaryDirectory(prefix='cfd-rust-exports-') as temp:
        destination = a.keep.resolve() if a.keep else Path(temp)
        for case in sorted((ROOT/'examples').glob('*.json')):
            out = destination/case.stem
            result = subprocess.run([str(binary), 'run', str(case), '--out', str(out)], check=True, capture_output=True, text=True)
            manifest = json.loads(result.stdout)
            assert manifest['execution_completed'] and not manifest['complete_clone']
            for filename in manifest['files']:
                path = out/filename
                assert path.is_file()
                if path.suffix == '.vtu':
                    validate_vtu(path, a.vtk_reader)
                    checked += 1
                elif path.suffix == '.pvd':
                    for item in ET.parse(path).findall('Collection/DataSet'):
                        assert math.isfinite(float(item.attrib['timestep']))
                        assert (path.parent/item.attrib['file']).is_file()
    runtime = 'not_executed'
    if a.vtk_reader:
        from vtkmodules.vtkCommonCore import vtkVersion
        runtime = 'VTK XML UnstructuredGrid reader ' + vtkVersion.GetVTKVersion()
    print(json.dumps({'case_examples': len(list((ROOT/'examples').glob('*.json'))), 'vtu_files_checked': checked,
                      'reader_runtime': runtime, 'paraview_application': 'not_executed'}))

if __name__ == '__main__':
    main()
