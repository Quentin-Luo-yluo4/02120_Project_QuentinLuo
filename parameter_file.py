"""Read model parameters from a text file, one per line:

    r = 1
    gamma_m = 0.5
    beta = 10
    gamma_p = 0.05
    t_max = 300

    # optional: the protein represses its own gene
    feedback = 1
    K = 200
    n = 2
    r_basal = 0

Rules:
  - each line is  name = value ; spaces around '=' are optional
  - blank lines are skipped, and '#' starts a comment (whole line or end of line)
  - names are the ones in PARAMETERS below (case-sensitive); each may appear once
  - parameters left out are simply not returned (the caller keeps its own value)
  - feedback = 1 turns feedback on, feedback = 0 turns it off; without that line,
    giving K, n or r_basal turns it on, and otherwise the setting is left as it was
Values are only checked to be ordinary numbers here; whether they make sense
(e.g. r >= 0) is checked by the model itself when it runs.

Use from Python:   parse_parameters(open("params.txt").read())
                   or read_parameter_file("params.txt")
"""

import math

PARAMETERS = {
    "r": "transcription rate",
    "gamma_m": "mRNA degradation rate",
    "beta": "translation rate",
    "gamma_p": "protein degradation rate",
    "t_max": "how long to simulate",
    "K": "feedback: protein level giving half repression",
    "n": "feedback: Hill coefficient",
    "r_basal": "feedback: rate at full repression",
    "feedback": "1 = protein represses its own gene, 0 = no feedback",
}
FEEDBACK_PARAMETERS = ("K", "n", "r_basal")


def parse_parameters(text):
    """Turn the file's text into a dict {name: float}.

    Raises ValueError naming the line number for anything it can't read.
    """
    values = {}
    for number, line in enumerate(text.splitlines(), start=1):
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        if "=" not in line:
            raise ValueError(f"line {number}: expected 'name = value', got {line!r}")
        name, value_text = (part.strip() for part in line.split("=", 1))
        if name not in PARAMETERS:
            raise ValueError(f"line {number}: unknown parameter {name!r}; "
                             f"allowed names are {', '.join(PARAMETERS)}")
        if name in values:
            raise ValueError(f"line {number}: {name} is given more than once")
        try:
            value = float(value_text)
        except ValueError:
            raise ValueError(f"line {number}: value of {name} is not a number: "
                             f"{value_text!r}") from None
        if not math.isfinite(value):
            raise ValueError(f"line {number}: value of {name} must be a finite number")
        if name == "feedback" and value not in (0.0, 1.0):
            raise ValueError(f"line {number}: feedback must be 0 (off) or 1 (on), got {value_text}")
        values[name] = value
    if not values:
        raise ValueError("no parameters found in the file")
    return values


def uses_feedback(values):
    """True/False if the file turns feedback on/off, None if it doesn't say.

    An explicit 'feedback = 0/1' wins; otherwise giving K, n or r_basal means on.
    """
    if "feedback" in values:
        return values["feedback"] == 1.0
    if any(name in values for name in FEEDBACK_PARAMETERS):
        return True
    return None


def read_parameter_file(path):
    with open(path, encoding="utf-8") as f:
        return parse_parameters(f.read())


if __name__ == "__main__":
    import sys
    print(read_parameter_file(sys.argv[1]))
