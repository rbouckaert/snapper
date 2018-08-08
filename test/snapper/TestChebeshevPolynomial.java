package snapper;

import java.util.Arrays;

import org.junit.Test;

import junit.framework.TestCase;

public class TestChebeshevPolynomial extends TestCase {
	
	@Test
	public void testChebeshevPolynomial() {
		int N = 33;
		ChebyshevPolynomial p = new ChebyshevPolynomial(N);
		p.init(60, 100);
		
		// store function values
		double [] f = new double[N];
		System.arraycopy(p.f, 0, f, 0, f.length);
		
		// convert factors to function values
		p.aToF();
		
		// check for largest error in function
		double maxErr = 0;
		int iErr = -1;
		for (int i = 0; i < f.length; i++) {
			if (Math.abs(f[i] - p.f[i]) > maxErr) {
				maxErr = Math.abs(f[i] - p.f[i]);
				iErr = i;
			}
		}
		System.out.println("Max error: " + maxErr + "\nat entry " + iErr);
		assertTrue(Math.abs(maxErr) < 1e-14);
	}

	
	@Test
	public void testChebeshevPolynomial2() {

		int N = 33;
		ChebyshevPolynomial p = new ChebyshevPolynomial(N);
		p.init(10, 32);

		// store function values
		double [] f = new double[N];
		System.arraycopy(p.f, 0, f, 0, f.length);
		System.out.println(Arrays.toString(f));
		double [] fin = new double[] {				
				0.00000000e+00,4.00424050e-19,3.49498706e-13,8.55254801e-10
				,1.75230607e-07,8.70224070e-06,1.67902916e-04,1.61580650e-03
				,8.98484427e-03,3.16714850e-02,7.52116961e-02,1.25324508e-01
				,1.50543597e-01,1.32619063e-01,8.64907774e-02,4.18856963e-02
				,1.50204264e-02,3.95490195e-03,7.53575618e-04,1.01716090e-04
				,9.44244301e-06,5.79210760e-07,2.22457764e-08,4.97232577e-10
				,5.84647512e-12,3.13273808e-14,6.18871196e-17,3.23947049e-20
				,2.57924641e-24,1.10157019e-29,2.42664019e-37,1.56377455e-50
				,0.00000000e+00				
				};
		for (int i = 0; i < N; i++) {
			System.out.println(fin[i] - f[i]);
		}
		double [] a = new double [] {
				2.10739037e-02,1.53264754e-02,-2.88498361e-02,-3.33576230e-02
				,1.80311475e-03,2.73147519e-02,1.59228114e-02,-9.11155969e-03
				,-1.58167632e-02,-3.70276402e-03,6.92228525e-03,5.84414927e-03
				,-2.73424845e-04,-2.87309624e-03,-1.38012637e-03,4.56719059e-04
				,7.45906917e-04,2.06888613e-04,-1.39383680e-04,-1.26441912e-04
				,-2.08314301e-05,2.15809534e-05,1.42207415e-05,1.88922381e-06
				,-1.72975078e-06,-1.02661462e-06,-1.86219733e-07,4.60513342e-08
				,3.70425070e-08,1.07994016e-08,1.79057434e-09,1.67866354e-10
				,6.99444148e-12	
		};
		System.out.println(Arrays.toString(p.a));
		for (int i = 0; i < N; i++) {
			System.out.println(a[i] - p.a[i]);
		}
		
		p.aToF();
		System.out.println(Arrays.toString(f));
		double [] fout = new double[]{
				
				-3.49721033e-12,3.49721422e-12,-3.14772324e-12,8.58752026e-10
				,1.75227109e-07,8.70224420e-06,1.67902913e-04,1.61580650e-03
				,8.98484426e-03,3.16714850e-02,7.52116961e-02,1.25324508e-01
				,1.50543597e-01,1.32619063e-01,8.64907774e-02,4.18856963e-02
				,1.50204264e-02,3.95490196e-03,7.53575614e-04,1.01716093e-04
				,9.44243951e-06,5.79214258e-07,2.22422792e-08,5.00729789e-10
				,2.34926488e-12,3.52854773e-12,-3.49715959e-12,3.49722334e-12
				,-3.49722106e-12,3.49723254e-12,-3.49722524e-12,3.49722640e-12
				,-3.49721727e-12
		};
		
		for (int i = 0; i < N; i++) {
			System.out.println(fout[i] - p.f[i]);
			assertEquals(fout[i], p.f[i], 5e-10);
		}
		
	}
}
