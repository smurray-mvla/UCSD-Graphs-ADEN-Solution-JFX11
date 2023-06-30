package week3example;



public class Coordinate {
	private int row;
	private int col;
	
	public Coordinate() {
		this.row = 0;
		this.col = 0;
	}
	
	public Coordinate(int row, int col) {
		this.row = row;
		this.col = col;
	}

	public int getRow() {
		return row;
	}
	
	public int getCol() {
		return col;
	}
	
}
