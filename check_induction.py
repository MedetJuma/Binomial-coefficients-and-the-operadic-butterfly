file = open("Xplus_lambda_GroebnerBasis.txt")

content = file.read().split("\n")

leading_terms = []
rest = []

for line in content:
    leading_terms.append(line.split("->")[0])

for i in range(len(leading_terms)):
    term = leading_terms[i]
    if not ("lambda" in term):
        continue
    # Check if term contains a Groebner basis leading term as fragment
    if "t(* x(* y" in term: # contains (33)
        continue
    if "lambda(* y" in term: # contains (26)
        continue
    if "lambda(* z" in term: # contains (15)
        continue
    if "x(lambda" in term: # contains (17)
        continue
    if "x(* y(* lambda(* t" in term: # reducible by x(* y(* lambda(* t(* *)))) (the inductive hypothesis)
        continue
    if "x(* y(* x(* lambda(* t" in term: # reducible by x(* y(* x(* lambda(* t(* *)) (the inductive hypothesis)
        continue

    print(content[i])

# The rules that were not reduced:
# x(* y(* lambda(* x(* t(* *)))))  ->  x(* t(* lambda(* x(* t(* *)))))
# x(* y(* x(* lambda(* x(* t(* *))))))  ->  x(* t(* x(* lambda(* x(* t(* *))))))
# x(* y(* lambda(* x(* x(* t(* *))))))  ->  x(* t(* lambda(* x(* x(* t(* *))))))
# x(* y(* lambda(* x(* x(* x(* t(* *)))))))  ->  x(* t(* lambda(* x(* x(* x(* t(* *)))))))
# x(* y(* x(* lambda(* x(* x(* t(* *)))))))  ->  x(* t(* x(* lambda(* x(* x(* t(* *)))))))
# x(* y(* x(* x(* lambda(* x(* t(* *)))))))  ->  x(* t(* x(* x(* lambda(* x(* t(* *)))))))
# x(* y(* lambda(* x(* x(* x(* x(* t(* *))))))))  ->  x(* t(* lambda(* x(* x(* x(* x(* t(* *))))))))
# x(* y(* x(* x(* x(* lambda(* x(* t(* *))))))))  ->  x(* t(* x(* x(* x(* lambda(* x(* t(* *))))))))
# x(* y(* x(* x(* lambda(* x(* x(* t(* *))))))))  ->  x(* t(* x(* x(* lambda(* x(* x(* t(* *))))))))
# x(* y(* x(* x(* lambda(* t(* x(* t(* *))))))))  ->  x(* t(* x(* x(* lambda(* t(* x(* t(* *))))))))
# x(* y(* x(* lambda(* x(* x(* x(* t(* *))))))))  ->  x(* t(* x(* lambda(* x(* x(* x(* t(* *))))))))
#
# These rules are precisely the elements R_p for various p.
