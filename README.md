This is the addendum to the paper "Binomial coefficients and the operadic butterfly."

The folders [Xminus](./Xminus/) and [Xplus](./Xplus/) contain the Haskell theory files "input.txt" for the operads $\mathcal X^-$ and $\mathcal X^+$, respectively. One can execute the Gröbner basis computation using the [Haskell operadic Gröbner bases calculator](https://irma.math.unistra.fr/~dotsenko/Operads.html). The file [Xminus_GroebnerBasis.txt](./Xminus/Xminus_GroebnerBasis.txt) contains the fully reduced Gröbner basis of $\mathcal X^-$, and [Xplus_GroebnerBasis.txt](./Xplus/Xplus_GroebnerBasis.txt) contains the Gröbner basis elements of $\mathcal X^+$ up to arity $6$.

Consider an operad $\mathcal X^+$ with an adjoined binary operation $\lambda$. The induction verification procedure files are contained in [Verification_Xplus](./Verification_Xplus/):

- The "input.txt" file is the Haskell theory file,
- The [Xplus_lambda_GroebnerBasis.txt](./Verification_Xplus/Xplus_lambda_GroebnerBasis.txt) file contains the Gröbner basis of the operad up to arity $9$,
- The [check_induction.py](./Verification_Xplus/check_induction.py) file contains the Python script that verifies that the rules arising from $$x(a_1, y(a_2, \lambda(a_3, x(a_4, t(a_5, a_6))))) \mapsto x(a_1, t(a_2, \lambda(a_3, x(a_4, t(a_5, a_6)))))$$ are reduced,, except $R_p$ for various $p$.
