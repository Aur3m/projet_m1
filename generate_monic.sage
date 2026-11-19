# géné polynomes moniques pr darboux (sage)
# TODO: nettoyer un peu si y a le tps


def generate_monomials(n_vars, degree):
    # genere les multidegres, recursif
    if n_vars == 1:
        return [(degree,)]

    monomials = []
    for first_power in range(degree + 1):
        remaining_degree = degree - first_power
        remaining_monomials = generate_monomials(n_vars - 1, remaining_degree)
        for rest in remaining_monomials:
            monomials.append((first_power,) + rest)

    return monomials


#fonction plus utilisée
def monomial_to_string(multideg, vars_list):
    # converti en string pour l'affichage
    terms = []
    for var, power in zip(vars_list, multideg):
        if power > 0:
            if power == 1:
                terms.append(str(var))
            else:
                terms.append(str(var) + "^" + str(power))

    if not terms:
        return "1"
    return "*".join(terms)


def generate_monic_polynomials_sage(var_names, max_degree, coeff_prefix='c'):
    # gros bloc pr fabriquer les polys moniques avec les c_i indéterminés
    n_vars = len(var_names)

    # 1. On liste d'abord tous les noms de coefficients nécessaires au total
    temp_coeff_names = []
    for target_degree in range(1, max_degree + 1):
        leading_monomials = generate_monomials(n_vars, target_degree)
        for _ in leading_monomials:
            for lower_deg in range(target_degree):
                lower_monomials = generate_monomials(n_vars, lower_deg)
                for _ in lower_monomials:
                    temp_coeff_names.append(f"{coeff_prefix}_{len(temp_coeff_names)}")

    all_vars = var_names.copy() + temp_coeff_names
    R = PolynomialRing(QQ, all_vars, order='lex')

    sage_vars = R.gens()
    var_dict = {var_names[i]: sage_vars[i] for i in range(n_vars)}
    
    # Dictionnaire associant le nom du coefficient à sa variable dans l'anneau R
    coeff_dict = {temp_coeff_names[i]: sage_vars[n_vars + i] for i in range(len(temp_coeff_names))}

    all_monomials = set()
    all_coefficients = set(temp_coeff_names)
    polynomials = []
    coeff_global_index = 0

    for deg in range(max_degree + 1):
        for monomial in generate_monomials(n_vars, deg):
            all_monomials.add(monomial)

    for target_degree in range(1, max_degree + 1):
        leading_monomials = generate_monomials(n_vars, target_degree)

        for leading_monomial in leading_monomials:
            poly_terms = []

            # leading term a coeff 1
            leading_term_sage = R(1)
            for i, power in enumerate(leading_monomial):
                for _ in range(power):
                    leading_term_sage *= sage_vars[i]

            poly_terms.append({
                'monomial': leading_monomial,
                'coefficient': 1,
                'coeff_var': None,
                'sage_term': leading_term_sage
            })

            # le reste avec des coeff indéterminés c_i
            for lower_deg in range(target_degree):
                lower_monomials = generate_monomials(n_vars, lower_deg)
                for monomial in lower_monomials:
                    coeff_var_name = temp_coeff_names[coeff_global_index]
                    coeff_var = coeff_dict[coeff_var_name]
                    coeff_global_index += 1

                    mono_sage = R(1)
                    for i, power in enumerate(monomial):
                        for _ in range(power):
                            mono_sage *= sage_vars[i]

                    poly_terms.append({
                        'monomial': monomial,
                        'coefficient': coeff_var,
                        'coeff_var': coeff_var_name,
                        'sage_term': mono_sage
                    })

            poly_sage = sum(term['coefficient'] * term['sage_term'] for term in poly_terms)

            polynomials.append({
                'degree': target_degree,
                'leading_monomial': leading_monomial,
                'terms': poly_terms,
                'undetermined_coefficients': [t['coeff_var'] for t in poly_terms if t['coeff_var']],
                'sage_poly': poly_sage
            })

    return {
        'ring': R,
        'var_dict': var_dict,
        'coeff_dict': coeff_dict,
        'monomials': sorted(list(all_monomials)),
        'polynomials': polynomials,
        'coefficients': sorted(list(all_coefficients)),
        'vars': var_names
    }


if __name__ == "__main__":
    # test pr verif si ca crash pas
    vars_list = ['x', 'y', 'z']
    max_deg = 3

    print("="*50)
    print("TEST GENERATION POLY MONIQUE")
    print("="*50)

    result = generate_monic_polynomials_sage(vars_list, max_deg)

    polynomials = result['polynomials']
    coefficients = result['coefficients']
    R = result['ring']

    print(f"nb polys générés : {len(polynomials)}")
    print(f"nb coeffs : {len(coefficients)}")

    for i, poly in enumerate(polynomials[:19], 1):
        print(f"\n{i}. Degré {poly['degree']}")
        print(f"   Poly : {poly['sage_poly']}")