# methode extactique, premier proto
# ref : cheze-combot [1] + christopher-llibre-pereira 2007 pour la propriete de divisibilite
# champ X = P d/dx + Q d/dy avec P, Q dans QQ[x,y]
# W = <m_1,...,m_N> (monomes de deg <= n), E_W = det(X^j(m_i)), j = 0..N-1
# si f darboux et f dans W alors f | E_W
# attention si E_W = 0 : il y a une integrale premiere rationnelle de deg <= n, et la ca sert a rien

R.<x,y> = PolynomialRing(QQ)

# monomes de deg total <= n, par deg croissant
def base_monomes(n):
    return [x^i * y^(k - i) for k in range(n + 1) for i in range(k + 1)]

# derivee le long du champ
def derivation(P, Q):
    return lambda f: P * f.derivative(x) + Q * f.derivative(y)

# ligne j = X^j applique a chaque elt de W -> matrice dim W x dim W
def matrice_extactique(P, Q, W):
    X = derivation(P, Q)
    lignes = [list(W)]
    for _ in range(1, len(W)):
        lignes.append([X(f) for f in lignes[-1]])
    return matrix(R, lignes)

# E_W, soit avec W donne soit avec tous les monomes de deg <= n
def polynome_extactique(P, Q, n=None, W=None):
    if W is None:
        W = base_monomes(n)
    return matrice_extactique(P, Q, W).determinant()

# on factorise E_n et on garde les facteurs qui sont vraiment darboux
# (tous les facteurs de E_n le sont pas !! cf lotka-volterra n=1)
# renvoie (E_n, [(f, cofacteur), ...]) ou (0, None) si E_n = 0
def candidats_darboux(P, Q, n):
    E = polynome_extactique(P, Q, n)
    if E == 0:
        return E, None
    X = derivation(P, Q)
    res = []
    for f, _ in E.factor():
        if f.degree() > n:
            continue
        g, r = X(f).quo_rem(f)   # darboux ssi le reste est nul
        if r == 0:
            res.append((f, g))
    return E, res

# pour voir combien ca coute quand n augmente (c'est le det qui explose)
def mesure_cout(P, Q, nmax):
    for n in range(1, nmax + 1):
        W = base_monomes(n)
        M = matrice_extactique(P, Q, W)
        t = cputime()
        E = M.determinant()
        t = cputime(t)
        deg = E.degree() if E != 0 else -1
        print("n = %d  taille %d x %d  temps %.2f s  deg E_n = %s"
              % (n, len(W), len(W), t, deg if deg >= 0 else "E_n = 0"))

# tests sur des cas ou on connait la reponse
if __name__ == "__main__":
    exemples = [
        ("rotation", -y, x),                      # x^2+y^2, cofacteur 0
        ("lotka-volterra", x*(1 - y), y*(x - 1)),  # x et y
    ]
    for nom, P, Q in exemples:
        print("==", nom)
        for n in [1, 2]:
            E, D = candidats_darboux(P, Q, n)
            if D is None:
                print("n = %d : E_n = 0 (integrale premiere rationnelle de deg <= %d)" % (n, n))
            else:
                print("n = %d : darboux trouves :" % n, D)
        mesure_cout(P, Q, 3)   # n=3 risque d'etre long
