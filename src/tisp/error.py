class TISPError(Exception):
    pass

class NotOptimalError(TISPError):
    pass

class UndefinedVariableError(TISPError):
    pass


class InfeasibleRelaxation(NotOptimalError):
    pass