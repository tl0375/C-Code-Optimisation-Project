#ifndef VTK_H
#define VTK_H

#ifdef __cplusplus
extern "C" {
#endif

void set_default_base();
void set_basename(char *base);
char *get_basename();
int write_checkpoint(int iteration);
int write_result();
int write_vtk(char* filename, int iters, double t);

#ifdef __cplusplus
}
#endif

#endif
