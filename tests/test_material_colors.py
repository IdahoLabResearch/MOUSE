import ast
from pathlib import Path

from matplotlib.colors import to_hex


ROOT = Path(__file__).parents[1]


def _literal_dict_assigned_to(path, variable_name):
    tree = ast.parse(path.read_text())
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        if any(isinstance(target, ast.Name) and target.id == variable_name for target in node.targets):
            return ast.literal_eval(node.value)
    raise AssertionError(f"Could not find {variable_name} in {path}")


def _material_database_keys():
    path = ROOT / 'core_design' / 'openmc_materials_database.py'
    tree = ast.parse(path.read_text())
    keys = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
            continue
        if node.func.attr != 'update' or not node.args or not isinstance(node.args[0], ast.Dict):
            continue
        if not isinstance(node.func.value, ast.Name) or node.func.value.id != 'materials_database':
            continue
        keys.update(ast.literal_eval(key) for key in node.args[0].keys)
    return keys


def test_every_database_material_has_one_unique_color():
    colors = _literal_dict_assigned_to(ROOT / 'core_design' / 'utils.py', 'MATERIAL_COLORS')

    assert set(colors) == _material_database_keys()
    resolved_colors = [to_hex(color) for color in colors.values()]
    assert len(resolved_colors) == len(set(resolved_colors))
