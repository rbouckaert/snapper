package snapper;

import java.util.Arrays;

import org.junit.Test;

import junit.framework.TestCase;

public class TestMatrixExponentiator extends TestCase {

	
	@Test
	public void testMatrixExponentiator() throws Exception {
		
		int N = 17;
		double [] a = new double[]{0.00661863, 0.00251134, -0.01136872, -0.00703671, 0.00858269, 0.00958962, -0.00366219, -0.01001947, -0.00015637, 0.00842034, 0.00338427, -0.0056598, -0.00479533, 0.0029283, 0.00431152, -0.00029069, -0.00434428};
		QMatrix Q = new QMatrix(N);
		Q.setQ(new double[]{.4,-0.5,0,0}, new double[]{0,0.5,-0.5,0});
		double [] b = MatrixExponentiator.expmvEuler(0.1, Q.Q, a, 0.0001);
		System.out.println(Arrays.toString(b));
		double [] b2 = new double[] {6.68075519e-03, 2.27250903e-03, -1.08318104e-02, -4.93046031e-03, 6.30544984e-03, 4.69595373e-03, -2.14186680e-03, -2.75532462e-03, 3.41168332e-04, 1.09133472e-03, 1.18947627e-04, -3.08169889e-04, -8.09840520e-05, 5.13696908e-05, 2.61197019e-05, -7.83636830e-06, -4.83855306e-06};
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
		
		int N = 17;
		double [] a = new double[]{0.00661863, 0.00251134, -0.01136872, -0.00703671, 0.00858269, 0.00958962, -0.00366219, -0.01001947, -0.00015637, 0.00842034, 0.00338427, -0.0056598, -0.00479533, 0.0029283, 0.00431152, -0.00029069, -0.00434428};
		QMatrix Q = new QMatrix(N);
		Q.setQ(new double[]{.4,-0.5,0,0}, new double[]{0,0.5,-0.5,0});
		double [] b = MatrixExponentiator.expmvRK4(0.1, Q.Q, a, 0.0001);
		System.out.println(Arrays.toString(b));
		double [] b2 = new double[] {6.68075519e-03, 2.27250903e-03, -1.08318104e-02, -4.93046031e-03, 6.30544984e-03, 4.69595373e-03, -2.14186680e-03, -2.75532462e-03, 3.41168332e-04, 1.09133472e-03, 1.18947627e-04, -3.08169889e-04, -8.09840520e-05, 5.13696908e-05, 2.61197019e-05, -7.83636830e-06, -4.83855306e-06};
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
		double [] a = new double[]{0.0, 0.00661863, 0.00251134, -0.01136872, -0.00703671, 0.00858269, 0.00958962, -0.00366219, -0.01001947, -0.00015637, 0.00842034, 0.00338427, -0.0056598, -0.00479533, 0.0029283, 0.00431152, -0.00029069, -0.00434428};
		QMatrix Q = new QMatrix(N);
		Q.setQ(new double[]{.4,-0.5,0,0}, new double[]{0,0.5,-0.5,0});
		double [] b = MatrixExponentiator.expmv(0.1, Q, a);
		System.out.println(Arrays.toString(b));
		double [] b2 = new double[] {0.0, 6.68075519e-03, 2.27250903e-03, -1.08318104e-02, -4.93046031e-03, 6.30544984e-03, 4.69595373e-03, -2.14186680e-03, -2.75532462e-03, 3.41168332e-04, 1.09133472e-03, 1.18947627e-04, -3.08169889e-04, -8.09840520e-05, 5.13696908e-05, 2.61197019e-05, -7.83636830e-06, -4.83855306e-06};
		// check for largest error in function
		double maxErr = 0;
		int iErr = -1;
		for (int i = 0; i < b.length; i++) {
			if (Math.abs(b[i] - b2[i]) > maxErr) {
				maxErr = Math.abs(b[i] - b2[i]);
				iErr = i;
			}
		}
		System.out.println("Max error M: " + maxErr + "\nat entry " + iErr);
	}
 }
