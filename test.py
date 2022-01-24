from distlink import COrbitData, MOID_fast, LC, test_peri_apo

o1 = COrbitData(7149.23810, 0.01951, 98.34212, 163.46357, 40.32555)
o2 = COrbitData(7254.82582, 0.02770, 98.10792, 523.30395-360, 10.76767)

dist = MOID_fast(o1, o2, 1e-10, 1e-10)
print(dist.distance, dist.distance_error)
print(dist.u1, dist.u1_error)
print(dist.u2, dist.u2_error)

link_coefs = LC(o1, o2, 0)
print(link_coefs.I, link_coefs.l, link_coefs.lmod, link_coefs.l2)

print(test_peri_apo(o1, o2, 100))