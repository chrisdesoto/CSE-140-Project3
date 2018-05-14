#define min(a,b) ((a < b) ? a : b)


// UNROLLING
void dgemm(int rows, int columns, float *A, float *C) {
  int unrolled = rows / 6;
  int loops;
  for (int col = 0; col < columns; ++col) {
    for(int row = 0; row < rows; ++row) {
      for(loops = 0; loops < unrolled; loops += 6) {
        C[loops + row * rows] += A[loops + col * rows] * A[row + col * rows];
        C[loops + 1 + row * rows] += A[loops + 1 + col * rows] * A[row + col * rows];
        C[loops + 2 + row * rows] += A[loops + 2 + col * rows] * A[row + col * rows];
        C[loops + 3 + row * rows] += A[loops + 3 + col * rows] * A[row + col * rows];
        C[loops + 4 + row * rows] += A[loops + 4 + col * rows] * A[row + col * rows];
        C[loops + 5 + row * rows] += A[loops + 5 + col * rows] * A[row + col * rows];
      }
      while (loops < rows) {
        C[loops + row * rows] += A[loops + rows * col] * A[row + rows * col];
      	loops++;
      }
    }
  }
}


/*
// TILING
void dgemm(int rows, int columns, float *A, float *C) {
	int tile_size = 32;
	for(int col = 0; col < columns; col += tile_size) {
		int tile_col_limit = min(col + tile_size, columns);
		for(int row = 0; row < rows; row += tile_size) {
			int tile_row_limit = min(row + tile_size, rows);
			for(int loop = 0; loop < rows; loop += tile_size) {
				int tile_loop_limit = min(loop + tile_size, rows);
				for(int tile_col = col; tile_col < tile_col_limit; ++tile_col) {
					for(int tile_row = row; tile_row < tile_row_limit; ++tile_row) {
						for(int tile_loop = loop; tile_loop < tile_loop_limit; ++tile_loop) {
							C[tile_row * rows + tile_loop] += A[rows * tile_col + tile_loop] * A[rows * tile_col + tile_row];
            }
          }
        }
			}
		}
	}
}
*/


// REORDERING
void dgemm(int m, int n, float *A, float *C) {
  for (int k = 0; k < n; k++) {
    for (int j = 0; j < m; j++) {
      for (int i = 0; i < m; i++) {
        C[i + j * m] += A[i + k * m] * A[j + k * m];
      }
    }
  }
}
