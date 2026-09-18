#!/usr/bin/env python3
"""
Nonsymmetric operadic Buchberger completion through arity 6.
Exact arithmetic (fractions.Fraction), deterministic, checkpointed.

It logs every nonzero S-polynomial remainder BEFORE monic normalization,
hence every divisor/denominator introduced by the GB algorithm.

Order: weighted path-deglex, w(t)=w(x)=1, w(y)=w(z)=4, t<y<x<z.
"""
from fractions import Fraction as F
from dataclasses import dataclass
from pathlib import Path
import pickle, json, hashlib, math

MAX_ARITY=6
OUT=Path("operad_buchberger_arity6_out"); OUT.mkdir(exist_ok=True)
OPS=("t","y","x","z"); W={'t':1,'x':1,'y':4,'z':4}; LR={'t':0,'y':1,'x':2,'z':3}

@dataclass(frozen=True)
class T:
    o:str|None; l:object=None; r:object=None
L=T(None)

def parse(s):
    s=''.join(s.split())
    def f(i):
        if s[i]=='*': return L,i+1
        o=s[i]; assert o in OPS and s[i+1]=='('
        a,j=f(i+2); b,j=f(j); assert s[j]==')'
        return T(o,a,b),j+1
    q,j=f(0); assert j==len(s); return q

def show(q): return '*' if q.o is None else f'{q.o}({show(q.l)} {show(q.r)})'
def ar(q): return 1 if q.o is None else ar(q.l)+ar(q.r)
def paths(q,w=()):
    return [w] if q.o is None else paths(q.l,w+(q.o,))+paths(q.r,w+(q.o,))
def total_weight(q):
    return 0 if q.o is None else W[q.o]+total_weight(q.l)+total_weight(q.r)
def wk(w): return (len(w),tuple(LR[a] for a in w))
def key(q): return (total_weight(q),tuple(wk(w) for w in paths(q)))
def lead(p): return max(p,key=key)
def clean(p): return {m:c for m,c in p.items() if c}
def add(p,q,a=F(1)):
    z=dict(p)
    for m,c in q.items():
        z[m]=z.get(m,F(0))+a*c
        if not z[m]: del z[m]
    return z

# Leaves of a pattern are distinct slots, left-to-right.
def match_at(pat,q):
    vals=[]
    def f(p,x):
        if p.o is None: vals.append(x); return True
        return x.o==p.o and f(p.l,x.l) and f(p.r,x.r)
    return vals if f(pat,q) else None

def inst(pat,vals):
    it=iter(vals)
    def f(p):
        if p.o is None:return next(it)
        return T(p.o,f(p.l),f(p.r))
    return f(pat)

def positions(q,pos=()):
    yield pos,q
    if q.o is not None:
        yield from positions(q.l,pos+(0,))
        yield from positions(q.r,pos+(1,))
def replace(q,pos,x):
    if not pos:return x
    return T(q.o,replace(q.l,pos[1:],x),q.r) if pos[0]==0 else T(q.o,q.l,replace(q.r,pos[1:],x))

def divisors(q,lm):
    ans=[]
    for pos,sub in positions(q):
        vals=match_at(lm,sub)
        if vals is not None: ans.append((pos,vals))
    return ans

def apply_rule_to_monomial(q,rule,pos,vals):
    # rule is monic lm + tail = 0, so lm -> -tail
    z={}
    for m,c in rule['poly'].items():
        if m==rule['lm']: continue
        n=replace(q,pos,inst(m,vals))
        z[n]=z.get(n,F(0))-c
    return clean(z)

def reductions(q,G):
    out=[]
    for gi,g in enumerate(G):
        for pos,vals in divisors(q,g['lm']):
            out.append((gi,pos,vals,apply_rule_to_monomial(q,g,pos,vals)))
    return out

def normal_form(p,G):
    p=dict(p); trace=[]
    while p:
        reducible=[]
        for m in p:
            rr=reductions(m,G)
            if rr: reducible.append((m,rr[0]))
        if not reducible: break
        m,(gi,pos,vals,rhs)=max(reducible,key=lambda x:key(x[0]))
        a=p.pop(m)
        p=add(p,rhs,a)
        trace.append((m,a,gi,pos))
    return p,trace

def leaf_paths(q,pos=(),out=None):
    if out is None: out=[]
    if q.o is None: out.append(pos)
    else:
        leaf_paths(q.l,pos+(0,),out); leaf_paths(q.r,pos+(1,),out)
    return out

def graft_at_leaf(q,i,x):
    ps=leaf_paths(q); return replace(q,ps[i],x)

def tree_positions(q,pos=()):
    yield pos
    if q.o is not None:
        yield from tree_positions(q.l,pos+(0,))
        yield from tree_positions(q.r,pos+(1,))

def subtree(q,pos):
    for a in pos:
        q=q.l if a==0 else q.r
    return q

def unify_patterns(a,b):
    """
    Least common refinement of two tree patterns whose leaves are wildcards.
    Returns the unique refinement when root positions are identified.
    """
    if a.o is None: return b
    if b.o is None: return a
    if a.o != b.o: return None
    l=unify_patterns(a.l,b.l)
    if l is None: return None
    r=unify_patterns(a.r,b.r)
    if r is None: return None
    return T(a.o,l,r)

def overlap_at_position(outer,inner,pos):
    """
    Identify the root of 'inner' with the vertex/leaf at pos in 'outer',
    unifying the two patterns there.  Leaves are slots and may be refined.
    """
    sub=subtree(outer,pos)
    u=unify_patterns(sub,inner)
    if u is None: return None
    return replace(outer,pos,u)

def internal_positions(q,pos=()):
    if q.o is None: return
    yield pos
    yield from internal_positions(q.l,pos+(0,))
    yield from internal_positions(q.r,pos+(1,))

def internal_vertex_set(q,pos=()):
    s=set()
    if q.o is not None:
        s.add(pos)
        s |= internal_vertex_set(q.l,pos+(0,))
        s |= internal_vertex_set(q.r,pos+(1,))
    return s

def occurrence_internal_vertices(pat,pos=()):
    """Positions of labelled/internal vertices of an occurrence of pat rooted at pos."""
    s=set()
    def f(q,p):
        if q.o is None:return
        s.add(p)
        f(q.l,p+(0,)); f(q.r,p+(1,))
    f(pat,pos); return s

def overlaps(lm1,lm2,maxarity):
    """
    Genuine minimal tree overlaps.

    Put the root of one leading monomial at every position of the other and
    unify the two tree patterns (leaves are wildcards).  Retain only common
    multiples in which the occurrences share an internal vertex and the union
    of their internal vertices is the whole common multiple.  This is the
    nonsymmetric tree analogue of a minimal common multiple.
    """
    ans=set()

    def try_oriented(a,b): # b rooted at pos of a
        for pos in tree_positions(a):
            q=overlap_at_position(a,b,pos)
            if q is None or ar(q)>maxarity: continue
            A=occurrence_internal_vertices(a,())
            B=occurrence_internal_vertices(b,pos)
            if not (A & B): continue
            if A | B != internal_vertex_set(q): continue
            # final divisibility sanity check
            if divisors(q,lm1) and divisors(q,lm2):
                ans.add(q)

    try_oriented(lm1,lm2)
    try_oriented(lm2,lm1)
    return sorted(ans,key=key)

def spolys_for_common_multiple(q,i,j,G):
    """
    Every occurrence of LM_i and LM_j in q gives two one-step reductions.
    Pair them; identical occurrence when i=j is skipped.
    """
    A=[]
    for pos,vals in divisors(q,G[i]['lm']):
        A.append((i,pos,vals,apply_rule_to_monomial(q,G[i],pos,vals)))
    B=[]
    for pos,vals in divisors(q,G[j]['lm']):
        B.append((j,pos,vals,apply_rule_to_monomial(q,G[j],pos,vals)))
    out=[]
    for a in A:
        for b in B:
            if i==j and a[1]==b[1]: continue
            s=add(a[3],b[3],F(-1))
            if s: out.append((a,b,s))
    return out

REL_BASE=[
[("x(x(* *) *)",1),("x(* x(* *))",-1),("x(* t(* *))",-1)],
[("x(x(* *) *)",1),("x(* y(* *))",-1),("x(* z(* *))",-1)],
[("x(z(* *) *)",1),("z(* x(* *))",-1),("z(* t(* *))",-1)],
[("z(x(* *) *)",1),("z(* z(* *))",-1),("z(* y(* *))",-1)],
[("z(z(* *) *)",1),("z(* z(* *))",-1),("z(* y(* *))",-1)],
[("x(t(* *) *)",1),("t(* x(* *))",-1)],
[("x(t(* *) *)",1),("t(* z(* *))",-1)],
[("x(y(* *) *)",1),("y(* x(* *))",-1)],
[("z(t(* *) *)",1),("y(* z(* *))",-1)],
[("z(y(* *) *)",1),("y(* z(* *))",-1)],
[("t(x(* *) *)",1),("t(t(* *) *)",1),("t(* t(* *))",-1)],
[("t(x(* *) *)",1),("t(t(* *) *)",1),("t(* y(* *))",-1)],
[("t(z(* *) *)",1),("t(y(* *) *)",1),("y(* t(* *))",-1)],
[("y(x(* *) *)",1),("y(t(* *) *)",1),("y(* y(* *))",-1)],
[("y(z(* *) *)",1),("y(y(* *) *)",1),("y(* y(* *))",-1)],
]

def relations_for(name):
    if name=="Xminus":
        r16=[("y(z(* *) *)",1),("y(x(* *) *)",-1),("x(* t(* *))",1),("x(* y(* *))",-1)]
    elif name=="Xplus":
        r16=[("y(z(* *) *)",1),("y(x(* *) *)",-1),("x(* t(* *))",-1),("x(* y(* *))",1)]
    else:
        raise ValueError(name)
    return REL_BASE+[r16]


def mkrel(a):
    z={}
    for s,c in a:
        m=parse(s); z[m]=z.get(m,F(0))+F(c)
    return clean(z)

def factors(n):
    n=abs(n); s=set(); d=2
    while d*d<=n:
        while n%d==0:s.add(d);n//=d
        d+=1
    if n>1:s.add(n)
    return s

def fmt(p):
    return ' + '.join(f'({c}) {show(m)}' for m,c in sorted(p.items(),key=lambda z:key(z[0]),reverse=True)) or '0'

def checkpoint(G,queue,done,events):
    with open(OUT/'checkpoint.pkl','wb') as f:
        pickle.dump(dict(G=G,queue=queue,done=done,events=events),f,pickle.HIGHEST_PROTOCOL)

def quadratic_rref(relations):
    rows=[mkrel(r) for r in relations]
    mons=sorted(set().union(*(r.keys() for r in rows)),key=key,reverse=True)
    R=[]; piv=[]
    for row in rows:
        row=dict(row)
        for pc,rr in zip(piv,R):
            if pc in row: row=add(row,rr,-row[pc])
        if not row: continue
        pc=next(m for m in mons if m in row)
        aa=row[pc]; row={m:c/aa for m,c in row.items()}
        for j,rr in enumerate(R):
            if pc in rr: R[j]=add(rr,row,-rr[pc])
        piv.append(pc); R.append(row)
    return sorted(zip(piv,R),key=lambda pr:key(pr[0]),reverse=True)

def run_operad(name):
    global OUT
    OUT=Path(f'operad_buchberger_{name}_out'); OUT.mkdir(exist_ok=True)
    pairs=quadratic_rref(relations_for(name))
    G=[dict(poly=row,lm=lm,arity=3,origin=f'quadratic RREF {k}')
       for k,(lm,row) in enumerate(pairs,1)]
    events=[]
    print(f'\n===== {name} =====')
    print('quadratic rank:',len(G))
    for i,g in enumerate(G): print(f'  G[{i}] LM {show(g["lm"])}')
    log=open(OUT/'events.txt','w'); stage_summary={}
    for n in range(4,MAX_ARITY+1):
        stage_primes=set(); passno=0; seen=set()
        while True:
            passno+=1; added=0; N=len(G); jobs=[]
            for i in range(N):
                for j in range(i,N):
                    for q in overlaps(G[i]['lm'],G[j]['lm'],n):
                        if ar(q)!=n: continue
                        for aa,bb,sp in spolys_for_common_multiple(q,i,j,G):
                            sig=(i,j,show(q),aa[1],bb[1])
                            if sig not in seen: jobs.append((i,j,q,aa,bb,sp,sig))
            jobs.sort(key=lambda J:(key(J[2]),J[0],J[1],J[3][1],J[4][1]))
            print(f'{name}: arity {n}, pass {passno}: {len(jobs)} unprocessed critical pairs; current rules={len(G)}',flush=True)
            for i,j,q,aa,bb,sp,sig in jobs:
                seen.add(sig); nf,tr=normal_form(sp,G)
                if not nf: continue
                lm=lead(nf); piv=nf[lm]
                if ar(lm)!=n: raise RuntimeError(f'arity changed {n}->{ar(lm)}')
                primes=sorted(factors(piv.numerator)|factors(piv.denominator)); stage_primes.update(primes)
                eno=len(events)+1
                ev=dict(kind='S',event=eno,pair=[i,j],common=show(q),arity=n,pivot=str(piv),
                        pivot_num=piv.numerator,pivot_den=piv.denominator,primes=primes,remainder=fmt(nf))
                events.append(ev); log.write(json.dumps(ev)+'\n'); log.flush()
                G.append(dict(poly={m:c/piv for m,c in nf.items()},lm=lm,arity=n,origin=f'S event {eno}'))
                added+=1; checkpoint(G,[],seen,events)
                print(f'  new G[{len(G)-1}] LM {show(lm)} pivot {piv} primes {primes}',flush=True)
            if added==0: break
        cnt=sum(g['arity']==n for g in G)
        stage_summary[str(n)]={'new_rules':cnt,'passes':passno,'pivot_primes':sorted(stage_primes)}
        print(f'{name}: ARITY {n} COMPLETE: {cnt} new rules; pivot primes {sorted(stage_primes)}')
    log.close()
    ps=set()
    for e in events: ps.update(e['primes'])
    summary={'operad':name,'rules':len(G),'rules_by_arity':{str(n):sum(g['arity']==n for g in G) for n in range(3,7)},
             'S_nonzero_events':len(events),'all_pivot_primes':sorted(ps),'stages':stage_summary}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    with open(OUT/'basis.txt','w') as f:
        for i,g in enumerate(G): f.write(f'G{i} [{g["origin"]}]: {fmt(g["poly"])}\n')
    checkpoint(G,[],set(),events); print(summary); return summary

def main():
    ans={name:run_operad(name) for name in ('Xminus','Xplus')}
    Path('operad_buchberger_Xminus_Xplus_summary.json').write_text(json.dumps(ans,indent=2)+'\n')
    print('\n===== COMBINED SUMMARY ====='); print(json.dumps(ans,indent=2))

if __name__=='__main__': main()
