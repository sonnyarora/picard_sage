from sage.rings.finite_rings import *
from sage.rings.all import PolynomialRing
from sage.rings.ideal import Ideal
from sage.rings.ideal import Ideal_generic
from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence
	
		

class PicardCurve():
	def __init__(self, f):
		self.F = f

	def base_field(self):
		P = self.F.parent()
		R = P.base_ring()
		return R;

	def get_eqn(self):
		return self.F;

	def genus(self):
		return 3;

	def __str__(self):
		return str(self.get_eqn());

class PicardJacobian():
	def __init__(self, C):
		# Check if Picard
		self.C = C;

	def random_point_typical(self, ff_q):
		C = self.get_curve()
		R.<z> = PolynomialRing(ff_q);
		temp_var = (C.get_eqn()).variable(0);
		f = 0;
		for i in range(0,5):
			f += ((C.get_eqn()).monomial_coefficient(temp_var^i))*z^i;

		f = - f;
		while true:
			v = 0;
			for i in range(0,3):
				v += (ff_q.random_element())*z^i;
			divs = divisors(v^3 - f);
			for u in divs:
				if u.leading_coefficient() == 1 and u.degree() >= v.degree() and u.degree() <= 3:
					return [u,v];

	def random_point(self, ff_q):
		T = TermOrder('M(3,4,0,1)')
		C = self.get_curve()
		ff_p = C.base_field()
		R.<x,y> = PolynomialRing(ff_q, 2, order=T)

		u, v = self.random_point_typical(ff_q);

		u_temp = 0;
		v_temp = 0; 
		for i in range(0,u.degree()+1):
			u_temp += (u.coefficients())[i]*x^(i);

		for i in range(0,v.degree()+1):
			v_temp += (v.coefficients())[i]*x^(i);

		return PicardJacobian_point(R.ideal([u_temp, y - v_temp]), self, R);


	def get_curve(self):
		return self.C;

	def get_curve_eqn(self):
		return (self.get_curve()).get_eqn();
		


	def reduce_point(I, F, R):
		B = I.groebner_basis();
		f = min_non_zero_fctn(R.ideal(B), R.ideal(F));
		J = R.ideal(f) + R.ideal(F);
		J = R.ideal(J.groebner_basis());
		J = J.quotient(I) + R.ideal(F);
		B = J.groebner_basis();
		return R.ideal(B);

	def min_non_zero_fctn(I,P):
		B = I.groebner_basis()
		i = len(B)

		while i > 0:
			f = B[i-1]
			if f.reduce(P) != 0:
				return f;
			i = i -1

		return 0;

# Based off Arita's code
class PicardJacobian_point():
	# Assumes pt is already Groebner reduced
	def __init__(self, pt, jac, poly_ring):
		self.pt = pt;
		self.jac = jac;
		self.poly_ring = poly_ring

	def get_jacobian(self):
		return self.jac;

	def get_ideal(self):
		return self.pt;

	def get_poly_ring(self):
		return self.poly_ring


	def is_zero(self):
		P = self.get_ideal();
		P_gens = P.groebner_basis();
		if len(P_gens) == 1:
			return true;

		return false;

	def frob(self, q):
		Q = self.get_ideal();
		R = self.get_poly_ring();
		gens = Q.groebner_basis();
		temp_gens = [];
		for gen in gens:
			temp_gen = 0;
			for mon in gen.monomials():

				temp_gen += R((gen.monomial_coefficient(mon))^q*mon);

			temp_gens.append(temp_gen);


		return PicardJacobian_point(R.ideal(temp_gens), self.get_jacobian(), R);



	def __add__(self, Q):
		F = (self.get_jacobian()).get_curve_eqn();
		R = self.get_poly_ring();
		I = (self.get_ideal())*Q.get_ideal() + R.ideal(F);

		J = reduce_point(I, F, R);
		J = reduce_point(J, F, R);

		return PicardJacobian_point(J, self.get_jacobian(), R);

	def __sub__(self, Q):
		
		return self + (-Q);

	def __neg__(self):
		F = (self.get_jacobian()).get_curve_eqn();
		R = self.get_poly_ring();
		I = self.get_ideal() + R.ideal(F);
		J = reduce_point(I, F, R);
		
		return PicardJacobian_point(J, self.get_jacobian(), R);



	def __mul__(self, n):
		r = PicardJacobian_point((self.get_poly_ring()).ideal(1), self.get_jacobian(), self.get_poly_ring());
		e = self;
		i = n;
		while i > 0:
			if (i % 2) == 1:				
				r = r + e
				i = (i-1)/2
			else:
				i = i/2
			if i >= 0:
				e = e + e;
		return r;


	def __str__(self):
		return str(self.get_ideal());

	def __eq__(self, Q):
		P = self.pt;
		Q = Q.get_ideal();
		Q_gens = Q.groebner_basis();
		P_gens = P.groebner_basis();
		if len(Q_gens) == len(P_gens):
			for i in range(0, len(P_gens)-1):
				if (Q_gens[i] in P_gens) == false:
					return false; # not the same gens
			
			return true;

		return false; #num gens not equal		


			


