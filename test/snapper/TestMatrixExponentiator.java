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
	
 }
