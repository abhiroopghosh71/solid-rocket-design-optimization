from problems.rocket_propellant_design import rocket


problem_list = {'rocket': rocket.RocketProblem}


def get_problem(name):
    if name not in problem_list:
        return None

    return problem_list[name.lower()]
