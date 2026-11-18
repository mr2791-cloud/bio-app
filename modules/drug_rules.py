def ghose_filter(mw, logp, atoms, molar_refractivity):
    """
    Ghose Filter:
    MW: 160-480
    LogP: -0.4 to 5.6
    Atom count: 20-70
    Molar refractivity: 40-130
    """

    checks = {
        "Molecular Weight": 160 <= mw <= 480,
        "LogP": -0.4 <= logp <= 5.6,
        "Atom Count": 20 <= atoms <= 70,
        "Molar Refractivity": 40 <= molar_refractivity <= 130
    }

    violations = sum(not value for value in checks.values())

    return {
        "checks": checks,
        "violations": violations,
        "passed": violations == 0
    }


def egan_rule(logp, tpsa):
    """
    Egan Rule:
    WLOGP <= 5.88
    TPSA <= 131.6
    """

    checks = {
        "LogP": logp <= 5.88,
        "TPSA": tpsa <= 131.6
    }

    violations = sum(not value for value in checks.values())

    return {
        "checks": checks,
        "violations": violations,
        "passed": violations == 0
    }
