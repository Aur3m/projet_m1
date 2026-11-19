# TODO: opti un peu les boucles d'index

def generate_monomials(n_vars, degree):
    # exposants homogènes de degré `degree` ; pour n_vars=2 : (0,d), (1,d-1), ..., (d,0)
    if n_vars == 1:
        return [(degree,)]
    monomials = []
    for first_power in range(degree + 1):
        remaining_degree = degree - first_power
        for rest in generate_monomials(n_vars - 1, remaining_degree):
            monomials.append((first_power,) + rest)
    return monomials


def init_ring(var_names, max_c=50, max_b=50, coeff_f='c', coeff_g='b'):
    # setup l'anneau avec les variables et les coeffs indéterminés
    c_names = [f"{coeff_f}_{i}" for i in range(max_c)]
    b_names = [f"{coeff_g}_{i}" for i in range(max_b)]
    all_vars = var_names.copy() + c_names + b_names

    # ordre lex pour groebner (x, y en tête, n'interviennent pas dans les équations)
    R = PolynomialRing(QQ, all_vars, order='lex')
    return R, c_names, b_names


def xy_term(sage_vars, m):
    term = 1
    for i, p in enumerate(m):
        term *= sage_vars[i]**p
    return term


def build_fi(R, sage_vars, c_names, k, leading_monomial):
    # f_i unitaire de terme dominant leading_monomial pour deg lex (x > y) :
    # tous les monômes de degré < k, et ceux de degré k strictement sous lm
    n_vars = 2
    lm = tuple(leading_monomial)
    f_i = R(xy_term(sage_vars, lm))

    c_index = 0
    for deg in range(k + 1):
        for m in generate_monomials(n_vars, deg):
            if deg == k and m >= lm:
                continue
            if c_index >= len(c_names):
                raise ValueError("max_c trop petit")
            f_i += R(c_names[c_index]) * xy_term(sage_vars, m)
            c_index += 1

    return f_i


def build_gi(R, sage_vars, b_names, n):
    # g_i générique de degré <= n
    n_vars = 2
    g_i = R(0)
    b_idx = 0
    for deg_g in range(n + 1):
        for m in generate_monomials(n_vars, deg_g):
            if b_idx >= len(b_names):
                raise ValueError("max_b trop petit")
            g_i += R(b_names[b_idx]) * xy_term(sage_vars, m)
            b_idx += 1
    return g_i


def get_xy_degree_and_lt(poly, x, y):
    # deg total en (x,y) et exposants du plus grand monôme deglex (x > y)
    # rq : dominant seulement génériquement, les coeffs dépendent des c_j
    xy_monos = {(m.degree(x), m.degree(y)) for m in poly.monomials()}
    if not xy_monos:
        return 0, 0, 0
    d = max(a + b for a, b in xy_monos)
    a, b = max(m for m in xy_monos if m[0] + m[1] == d)
    return d, a, b


#K[x,y] est hardcodé ici, la fonction n'est pas générique mais pas nécessaire pr l'instant
def extract_equations(E, R):
    # sort les eqns depuis E = f_i*g_i - D_fi : un coeff (polynôme en c, b) par monôme en x,y
    E_dict = {}
    for mono, coeff in E.dict().items():
        mono = tuple(mono)   # ETuple -> tuple
        xy_tuple = mono[:2]
        rest_tuple = mono[2:]
        coeff_poly = coeff * R({(0, 0) + rest_tuple: 1})
        E_dict[xy_tuple] = E_dict.get(xy_tuple, R(0)) + coeff_poly
    return E_dict


def split_gequal_fequal(E_dict, n, lt_fi_x, lt_fi_y):
    # separe entre gequal et fequal
    G_monos = set()
    for deg_g in range(n + 1):
        for m in generate_monomials(2, deg_g):
            G_monos.add((m[0] + lt_fi_x, m[1] + lt_fi_y))

    #genere tout les monômes qu'on obtient en multipliant LT(f_i) par un monôme quelconque de g_i

    gequal = []
    fequal = []
    for xy_tuple, eq in E_dict.items():
        if eq == 0:
            continue
        if xy_tuple in G_monos:
            # Si (x, y) ∈ G_monos :
            # Ce monôme vient du produit LT(f_i) · g_i
            gequal.append(eq)
        else:
            # Si (x, y) ∉ G_monos :
            # Ce monôme vient des termes mineurs de f_i ou de D(f_i)
            #C'est une équation qu'on résoudra avec Gröbner
            fequal.append(eq)

    return gequal, fequal


def ps1_groebner_modular(P, Q, var_names, max_degree):
    # main qui appele tout
    R, c_names, b_names = init_ring(var_names)
    sage_vars = R.gens()
    x, y = sage_vars[0], sage_vars[1]

    P_ring = R(P)
    Q_ring = R(Q)
    S_fg = []

    for k in range(1, max_degree + 1):
        print(f"--> Check degré k={k}")
        for lm in generate_monomials(2, k):
            lt_fi_x, lt_fi_y = lm[0], lm[1]

            f_i = build_fi(R, sage_vars, c_names, k, lm)

            # D(f_i)
            D_fi = P_ring * f_i.derivative(x) + Q_ring * f_i.derivative(y)
            if D_fi == 0:
                continue

            deg_Dfi, lt_Df_x, lt_Df_y = get_xy_degree_and_lt(D_fi, x, y)

            # check de divisibilité
            # Filtre 1 : LT(D(f_i)) doit être divisible par LT(f_i)
            if not (lt_fi_x <= lt_Df_x and lt_fi_y <= lt_Df_y):
                continue
            # Filtre 2 : deg(D(f_i)) doit être au moins k
            if deg_Dfi < k:
                continue

            # Calculer le degré de g_i
            n = deg_Dfi - k

            #construction de l'ansatz g_i
            g_i = build_gi(R, sage_vars, b_names, n)

            #extraction des equations
            E = f_i * g_i - D_fi
            E_dict = extract_equations(E, R)
            gequal, fequal = split_gequal_fequal(E_dict, n, lt_fi_x, lt_fi_y)

            # resolution par ideal/groebner
            eqs = gequal + fequal
            if not eqs:
                gb = []
            else:
                J = R.ideal(eqs)
                if R(1) in J:
                    continue  # système incompatible
                gb = J.groebner_basis()

            print(f"  -> Sol trouvée ! lt(fi) = x^{lt_fi_x} y^{lt_fi_y}")
            S_fg.append({
                'f_i': f_i,
                'g_i': g_i,
                'groebner': gb
            })

    return S_fg


if __name__ == "__main__":
    R_temp = PolynomialRing(QQ, ['x', 'y'])
    x, y = R_temp.gens()

    # (2x^2*y - x) dy/dx = 2xy^2 + y
    P = 2*x**2*y - x
    Q = 2*x*y**2 + y

    print("=" * 50)
    print("TEST PS1 MODULAIRE")
    print("=" * 50)

    resultats = ps1_groebner_modular(P, Q, ['x', 'y'], max_degree=2)

    print(f"\nNb solutions polynomiales (S_fg) : {len(resultats)}")
    for idx, sol in enumerate(resultats):
        print(f"\nSol #{idx+1}:")
        print(f"f_i générique : {sol['f_i']}")
        print(f"g_i générique : {sol['g_i']}")
        print(f"Base de Groebner : {sol['groebner']}")
