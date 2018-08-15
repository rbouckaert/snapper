package snapper;

import java.util.Arrays;

import org.junit.Test;

import junit.framework.TestCase;

public class TestMatrixExponentiator extends TestCase {
	int N = 17;
	double [] a = new double[]{0.00661863, 0.00251134, -0.01136872, -0.00703671, 0.00858269, 0.00958962, -0.00366219, -0.01001947, -0.00015637, 0.00842034, 0.00338427, -0.0056598, -0.00479533, 0.0029283, 0.00431152, -0.00029069, -0.00434428};
	double [] b2 = new double[] {6.68075519e-03, 2.27250903e-03, -1.08318104e-02, -4.93046031e-03, 6.30544984e-03, 4.69595373e-03, -2.14186680e-03, -2.75532462e-03, 3.41168332e-04, 1.09133472e-03, 1.18947627e-04, -3.08169889e-04, -8.09840520e-05, 5.13696908e-05, 2.61197019e-05, -7.83636830e-06, -4.83855306e-06};

	
	@Test
	public void testMatrixExponentiator() throws Exception {
		
		QMatrix Q = new QMatrix(N);
		Q.setQ(new double[]{.4,-0.5,0,0}, new double[]{0,0.5,-0.5,0});
		double [] b = MatrixExponentiator.expmvEuler(0.1, Q.Q, a, 0.001);
		System.out.println(Arrays.toString(b));
		// check for largest error in function
		double maxErr = 0;
		int iErr = -1;
		for (int i = 0; i < b.length; i++) {
			if (Math.abs(b[i] - b2[i]) > maxErr) {
				maxErr = Math.abs(b[i] - b2[i]);
				iErr = i;
			}
		}
		System.out.println("Max error Euler: " + maxErr + "\nat entry " + iErr);
	}
	
	@Test
	public void testMatrixExponentiator2() throws Exception {
		
		QMatrix Q = new QMatrix(N);
		Q.setQ(new double[]{.4,-0.5,0,0}, new double[]{0,0.5,-0.5,0});
		double [] b = MatrixExponentiator.expmvRK4(0.1, Q.Q, a, 0.001);
		System.out.println(Arrays.toString(b));
		// check for largest error in function
		double maxErr = 0;
		int iErr = -1;
		for (int i = 0; i < b.length; i++) {
			if (Math.abs(b[i] - b2[i]) > maxErr) {
				maxErr = Math.abs(b[i] - b2[i]);
				iErr = i;
			}
		}
		System.out.println("Max error RK4: " + maxErr + "\nat entry " + iErr);
	}
	
	@Test
	public void testMatrixExponentiator3() throws Exception {
		
		int N = 17;
		double [] a = new double[N+1];
		System.arraycopy(this.a, 0, a, 1, N);
		QMatrix Q = new QMatrix(N);
		Q.setQ(new double[]{.4,-0.5,0,0}, new double[]{0,0.5,-0.5,0});
		double [] b = MatrixExponentiator.expmv(0.1, Q, a);
		System.out.println(Arrays.toString(b));
		double [] b2 = new double[N+1];
		System.arraycopy(this.b2, 0, b2, 1, N);
		// check for largest error in function
		double maxErr = 0;
		int iErr = -1;
		for (int i = 0; i < b.length; i++) {
			if (Math.abs(b[i] - b2[i]) > maxErr) {
				maxErr = Math.abs(b[i] - b2[i]);
				iErr = i;
			}
		}
		System.out.println("Max error expokit: " + maxErr + "\nat entry " + iErr);
	}
	
	
	double [] a2 = new double[]{1.13078299e-04,-3.64183566e-05,-2.00360365e-04,9.22944904e-05
	       ,1.37472763e-04,-1.09621517e-04,-6.89952861e-05,9.19762926e-05
	       ,2.00984589e-05,-5.92555604e-05,2.96422515e-06,3.00218468e-05
	       ,-8.09975061e-06,-1.19199919e-05,5.92212653e-06,3.57768294e-06
	       ,-2.91045842e-06,-7.12537520e-07,1.08736266e-06,3.02266301e-08
	       ,-3.20207576e-07,4.40459701e-08,7.47387270e-08,-2.19525795e-08
	       ,-1.35550447e-08,6.56674018e-09,1.77754118e-09,-1.45889932e-09
	       ,-1.21575348e-10,2.54241133e-10,-1.31983020e-11,-3.15276178e-11
	       ,-1.57638089e-11};
	double [] a3 = new double[]{1.13078299592331E-4, -3.6418357807009555E-5, -2.0036036417602466E-4, 9.229449024964959E-5, 1.3747276267067577E-4, -1.0962151594729884E-4, -6.899528667937696E-5, 9.197629315280765E-5, 2.009845848857312E-5, -5.9255560178417786E-5, 2.9642252006610002E-6, 3.002184652850581E-5, -8.09975045367687E-6, -1.1919991931487868E-5, 5.9221265278692834E-6, 3.5776828690870717E-6, -2.9104583719774233E-6, -7.125374953547786E-7, 1.0873626300027483E-6, 3.02266360057486E-8, -3.202075720515647E-7, 4.404596933163639E-8, 7.473872819449082E-8, -2.195258017882337E-8, -1.3555046283621306E-8, 6.566741157131833E-9, 1.777542026880363E-9, -1.4589000954583843E-9, -1.2157552474646597E-10, 2.5424146516778113E-10, -1.319830837909165E-11, -3.152769802091402E-11, 5.6723277808211174E-12};
	                           
	double [] b3 = new double[]{1.22958683e-04,-4.22063593e-05,-1.90385477e-04,8.25912678e-05
	       ,8.73729113e-05,-5.91860621e-05,-2.22159841e-05,2.36905386e-05
	       ,2.11242554e-06,-5.81376962e-06,4.01804130e-07,9.01553186e-07
	       ,-1.59990993e-07,-8.81671434e-08,2.42726090e-08,5.19795264e-09
	       ,-2.13136228e-09,-1.52335365e-10,1.18207709e-10,-8.39580591e-13
	       ,-4.25567210e-12,2.50088857e-13,9.95655923e-14,-1.00094945e-14
	       ,-1.47545560e-15,2.16781038e-16,1.26405989e-17,-2.92491372e-18
	       ,-4.03500396e-20,2.51032811e-20,1.63718020e-22,-1.13509330e-22
	       ,-1.08997282e-23};
	
	@Test
	public void testMatrixExponentiator4() throws Exception {
		
		int N = 33;
		QMatrix Q = new QMatrix(N);
		Q.setQ(new double[]{.5,-1,0,0}, new double[]{0,.5,-.5,0});

//		// expokit uses arrays offset by 1, so reserve N+1 places
//		double [] a = new double[N+1];
//		System.arraycopy(this.a3	, 0, a, 1, N);
//		double [] b = MatrixExponentiator.expmv(0.1, Q, a);
//		System.out.println(Arrays.toString(b));
//		double [] b3 = new double[N+1];
//		System.arraycopy(this.b3, 0, b3, 1, N);
//		// check for largest error in function
//		double maxErr = 0;
//		int iErr = -1;
//		for (int i = 0; i < b.length; i++) {
//			if (Math.abs(b[i] - b3[i]) > maxErr) {
//				maxErr = Math.abs(b[i] - b3[i]);
//				iErr = i;
//			}
//		}
//		System.out.println("Max error expokit2: " + maxErr + "\nat entry " + iErr);

		double [] b = MatrixExponentiator.expmvRK4(0.1, Q.Q, a3, 0.001);
		System.out.println(Arrays.toString(b));
		// check for largest error in function
		double maxErr = 0;
		int iErr = -1;
		for (int i = 0; i < b.length; i++) {
			if (Math.abs(b[i] - b3[i]) > maxErr) {
				maxErr = Math.abs(b[i] - b3[i]);
				iErr = i;
			}
		}
		System.out.println("Max error RK4(2): " + maxErr + "\nat entry " + iErr);	}
	

	@Test
	public void testMatrixExponentiator5() {
		int N = 33;
		ChebyshevPolynomial p = new ChebyshevPolynomial(N);
		p.init(60, 100);

		// store function values
		double [] f = new double[N];
		System.arraycopy(p.f, 0, f, 0, f.length);
		System.out.println(Arrays.toString(f));
		double [] fin = new double[] {
				0.00000000e+000,9.81217458e-130,8.44740350e-094,5.52746404e-073
				,1.96058626e-058,2.23637075e-047,1.40561393e-038,2.20794099e-031
				,2.12897558e-025,2.24025170e-020,3.79277862e-016,1.35570782e-012
				,1.24328444e-009,3.37020772e-007,2.99255612e-005,9.36366872e-004
				,1.08438667e-002,4.78258544e-002,8.11173785e-002,5.24425537e-002
				,1.25564545e-002,1.05721928e-003,2.88859854e-005,2.27585665e-007
				,4.35710569e-010,1.57796090e-013,7.41815666e-018,2.50114523e-023
				,2.21775979e-030,7.81024821e-040,1.55163732e-053,2.18277324e-077
				,0.00000000e+000
				};
		for (int i = 0; i < N; i++) {
			System.out.println(fin[i] - f[i]);
		}
		double [] a = new double [] {
				6.46372099e-03,-2.55988950e-03,-1.16700845e-02,6.99254795e-03
				,8.30736893e-03,-9.62449434e-03,-3.89177042e-03,9.99997936e-03
				,-3.16413803e-04,-8.41597825e-03,3.32499151e-03,5.70139468e-03
				,-4.70808001e-03,-2.82391503e-03,4.62516343e-03,5.27632266e-04
				,-3.61099170e-03,8.51152558e-04,2.28218230e-03,-1.36076693e-03
				,-1.10635964e-03,1.27694765e-03,3.07749738e-04,-9.21548777e-04
				,1.02340235e-04,5.38834238e-04,-2.32248762e-04,-2.54010441e-04
				,2.18037296e-04,9.22966180e-05,-1.67916650e-04,-2.01820567e-05
				,-1.00910284e-05
		};
		System.out.println(Arrays.toString(p.a));
		for (int i = 0; i < N; i++) {
			System.out.println(a[i] + " " +  p.a[i] + " " + (a[i] - p.a[i]));
	//		assertEquals(a[i], p.a[i], 5e-10);
		}
		
		double [] exp = new double[]{
				7.01716030e-03,-2.94255425e-03,-1.11341912e-02,6.33565510e-03
				,5.31269743e-03,-5.42670684e-03,-1.12815346e-03,2.77378802e-03
				,-2.69896780e-04,-9.05679289e-04,2.87824983e-04,1.84696096e-04
				,-1.05203315e-04,-1.94253323e-05,2.28828711e-05,-5.28541986e-07
				,-3.22723820e-06,5.27386625e-07,2.93315541e-07,-8.94030807e-08
				,-1.53099813e-08,8.66637254e-09,1.51203823e-10,-5.45194196e-10
				,4.16486463e-11,2.28086803e-11,-3.75075989e-12,-6.22940046e-13
				,1.69970658e-13,1.01503240e-14,-2.73890476e-15,-7.26617451e-17
				,-6.97734076e-18
		};
		QMatrix Q = new QMatrix(N);
		Q.setQ(new double[]{.5,-1,0,0}, new double[] {0,0.5,-0.5,0});
		double [] b = MatrixExponentiator.expmvRK4(0.1, Q.Q, p.a, 0.001);
		for (int i = 0; i < N; i++) {
			System.out.println((b[i] - exp[i]));
			assertEquals(exp[i], b[i], 5e-8);
		}
	}

}
