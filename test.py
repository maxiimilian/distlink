from distlink import COrbitData, MOID_fast, LC, test_peri_apo

o1 = COrbitData(7149.23810, 0.01951, 1.7164, 2.853, 0.7038)
o2 = COrbitData(7254.82582, 0.02770, 1.7123, 2.85, 0.188)

dist = MOID_fast(o1, o2, 2e-15, 1e-15)
print(dist.distance, dist.distance_error)
print(dist.u1, dist.u1_error)
print(dist.u2, dist.u2_error)

link_coefs = LC(o1, o2, 0)
print(link_coefs.I, link_coefs.l, link_coefs.lmod, link_coefs.l2)

print(test_peri_apo(o1, o2, 100))