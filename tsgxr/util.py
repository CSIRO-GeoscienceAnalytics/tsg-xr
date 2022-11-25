from pathlib import Path


def rm_tree(pth: Path):
    """Recursively remove a folder."""
    if pth.exists():
        if pth.is_file():
            pth.unlink()
        else:
            for child in pth.iterdir():
                if child.is_file():
                    child.unlink()
                else:
                    rm_tree(child)
            pth.rmdir()
