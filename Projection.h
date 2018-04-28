struct Tile
{
	int fid; /* face id */
	int tid_h; /* Tile id in horizontal axis */
	int tid_w; /* Tile id in vertical axis */
};
int* getAdjacentTileCMP(int NO_TILE_H, int NO_TILE_W, int tid, int* N);
int* getAdjacentTileERP2(int NO_TILE_H, int NO_TILE_W, int tid, int* N);
int* getAdjacentTileERP(int NO_TILE_H, int NO_TILE_W, int tid, int* N);
int* getAdjacentTile(int NO_FACE, int NO_TILE_H, int NO_TILE_W, int tid, int *N);
int* extVmask(int* vmask,int NO_FACE, int NO_TILE_H, int NO_TILE_W, int ext_width);
void showVmask(int* vmask, int NO_FACE, int NO_TILE_H, int NO_TILE_W);
void showArray(int *arr, int N);
Tile* get_adj_tile_in_face(int NO_TILE_H, int NO_TILE_W, int tid_w, int tid_h, int* N);