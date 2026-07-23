"""Single source of truth for the genome assembly pipeline version."""

__version__ = "0.2.2"
__version_info__ = (0, 2, 2, None)


def parse_version(s):
    """Parse a PEP 440-lite version string.

    Accepted shape is MAJOR.MINOR.PATCH with optional aN, bN, or rcN suffix.
    """
    import re

    m = re.match(r"^(\d+)\.(\d+)\.(\d+)((?:a|b|rc)\d+)?$", s.strip())
    if m is None:
        raise ValueError(
            f"not a PEP 440 lite version: {s!r} "
            f"(expected MAJOR.MINOR.PATCH with optional aN/bN/rcN modifier)"
        )
    return (int(m.group(1)), int(m.group(2)), int(m.group(3)), m.group(4))


if __name__ == "__main__":
    print(__version__)
