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
		
		double [] M = new double[N*N];
		Q.getQ(new double[]{.4,-0.5,0,0}, new double[]{0,0.5,-0.5,0}, M);
		Q.prettyPrint(M);
		
	}
}
