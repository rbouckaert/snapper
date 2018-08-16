package snapper;

import java.util.Arrays;

import org.junit.Test;

import junit.framework.TestCase;

public class TestLikelihoodCore extends TestCase {

	@Test
	public void testRootIntegrator() {
		double [] F1 = new double[]{1.65068410e-04,-1.13050435e-04,-2.01721354e-04,1.89332689e-04
		          ,1.84667834e-05,-9.69752622e-05,3.24319878e-05,1.96353865e-05
		          ,-1.54266955e-05,1.79448453e-07,2.92513803e-06,-7.63365028e-07
		          ,-2.26080621e-07,1.35546193e-07,-3.76539278e-09,-1.05607964e-08
		          ,1.95391485e-09,3.40594276e-10,-1.44248374e-10,3.27552084e-12
		          ,4.85886819e-12,-6.07951773e-13,-6.87174125e-14,1.95762501e-14
		          ,-2.93384805e-16,-2.92277888e-16,2.41435326e-17,1.86244671e-18
		          ,-3.47137839e-19,2.26479866e-21,2.34658482e-21,-7.63583736e-23
		          ,-7.33230933e-24};
		double [] F2 = new double[]{1.22958683e-04,-4.22063593e-05,-1.90385477e-04,8.25912678e-05
		          ,8.73729113e-05,-5.91860621e-05,-2.22159841e-05,2.36905386e-05
		          ,2.11242554e-06,-5.81376962e-06,4.01804130e-07,9.01553186e-07
		          ,-1.59990993e-07,-8.81671434e-08,2.42726090e-08,5.19795264e-09
		          ,-2.13136228e-09,-1.52335365e-10,1.18207709e-10,-8.39580591e-13
		          ,-4.25567210e-12,2.50088857e-13,9.95655923e-14,-1.00094945e-14
		          ,-1.47545560e-15,2.16781038e-16,1.26405989e-17,-2.92491372e-18
		          ,-4.03500396e-20,2.51032811e-20,1.63718020e-22,-1.13509330e-22
		          ,-1.08997282e-23};
		int N = 33;
		ChebyshevPolynomial p1 = new ChebyshevPolynomial(N);
		System.arraycopy(F1, 0, p1.a, 0, N);
		p1.aToF();

		ChebyshevPolynomial p2 = new ChebyshevPolynomial(N);
		System.arraycopy(F2, 0, p2.a, 0, N);
		p2.aToF();
		
		double [] rootPartials = new double[N];
		for (int i = 0; i < N; i++) {
			rootPartials[i] = p1.f[i] * p2.f[i];
		}
		
		SnapperLikelihoodCore  core = new SnapperLikelihoodCore(3, 1, N);
		core.initialize(3, 1, 1, true, false);
		double [] outLogLikelihoods = new double[1];
		double [] frequencies = new double[N];
		Arrays.fill(frequencies, 1);
		core.calculateLogLikelihoods(rootPartials, frequencies, outLogLikelihoods);
		System.out.println("logL = " + outLogLikelihoods[0]);
		//    logL = -16.359225817728493 // trapezoid
		//    logL = -16.35761892714236  // Clenshaw-Curtis quadrature
		assertEquals(-16.357618928099296, outLogLikelihoods[0], 1e-4);
	
	}
}
