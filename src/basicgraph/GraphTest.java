package basicgraph;

import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

class GraphTest {

	@BeforeAll
	static void setUpBeforeClass() throws Exception {
	}

	@AfterAll
	static void tearDownAfterClass() throws Exception {
	}

	@BeforeEach
	void setUp() throws Exception {
	}

	@AfterEach
	void tearDown() throws Exception {
	}

	@Test
	void test() {
		GraphAdjList lst = new GraphAdjList();
		GraphAdjMatrix mat = new GraphAdjMatrix();
		
		for (int i = 0; i < 4; i++) {
			mat.addVertex();
			lst.addVertex();
		}
		
		mat.addEdge(0, 1);
		mat.addEdge(0, 2);
		mat.addEdge(1, 3);
		mat.addEdge(2, 1);
		mat.addEdge(2, 3);

		lst.addEdge(0, 1);
		lst.addEdge(0, 2);
		lst.addEdge(1, 3);
		lst.addEdge(2, 1);
		lst.addEdge(2, 3);

		for (int i = 0; i < 4; i++) {
			System.out.println("Matrix: Get Neighbors vertex = "+i+" returns "+mat.getNeighbors(i).toString());
			System.out.println("  List: Get Neighbors vertex = "+i+" returns "+lst.getNeighbors(i).toString());
		}
		for (int i = 0; i < 4; i++) {
			System.out.println("Matrix: Get Distance2 vertex = "+i+" returns "+mat.getDistance2(i).toString());
			System.out.println("  List: Get Distance2 vertex = "+i+" returns "+lst.getDistance2(i).toString());
		}

		assertTrue(true);
	}


}
