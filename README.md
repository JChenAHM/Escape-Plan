Escape-Plan
===========
Given a directed graph G = (V, E) where the vertices represent populated cities.

The vertices can be categorized in three classes:

i) X⊂V, X denotes the set of populated but ‘unsafe’ cities.

ii) S⊂V, S denotes the set of ‘safe cities’.

iii) N⊂V, N denotes the set of cities which are not categorized ‘safe’ or ‘unsafe’. N can be empty.

X ∩ S = S ∩ N = N ∩ X = ϕ.

X ∪ S ∪ N = V

Our goal is to solve the following problem:

If there exists a path from each ‘unsafe’ city (x ϵ X) to a ‘safe’ city (s ϵ S) such that no two paths have any common edge,
then we say there is an ‘Escape Plan’. (They can have common vertices.) Does an escape plan exist? 
