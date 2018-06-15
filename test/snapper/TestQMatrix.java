package snapper;

import org.junit.Test;

import junit.framework.TestCase;

public class TestQMatrix extends TestCase {

	@Test
	public void testQMatrix() {
		int N = 10;
		QMatrix Q = new QMatrix(N);
		Q.prettyPrint(Q.D);
		Q.prettyPrint(Q.S2D2);
		
		Q.setQ(new double[]{.4,-0.5,0,0}, new double[]{0,0.5,-0.5,0});
		double [] M = Q.Q;
		Q.prettyPrint(M);
		
		double [] expectedQ = new double [] {0., 0.3, -0.5, 0.9, -1., 1.5, -1.5, 2.1, -2., 2.7, 0., -0.5, 1.2, -1.5, 2.4, -2.5, 3.6, -3.5, 4.8, -4.5, 0., 0., -1.5, 1.8, -2., 3., -3., 4.2, -4., 5.4, 0., 0., 0., -3., 2.4, -2.5, 3.6, -3.5, 4.8, -4.5, 0., 0., 0., 0., -5., 3., -3., 4.2, -4., 5.4, 0., 0., 0., 0., 0., -7.5, 3.6, -3.5, 4.8, -4.5, 0., 0., 0., 0., 0., 0., -10.5, 4.2, -4., 5.4, 0., 0., 0., 0., 0., 0., 0., -14., 4.8, -4.5, 0., 0., 0., 0., 0., 0., 0., 0., -18., 5.4, 0., 0., 0., 0., 0., 0., 0., 0., 0., -22.5};
		for (int i = 0; i < N*N; i++) {
			assertEquals(expectedQ[i], M[i], 1e-13);
		}
	}
}
