#include <LAGraph.h>
#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

int LAGraph_msf(GrB_Matrix *result, GrB_Matrix A, bool sanitize, char *msg);

typedef struct {
    GrB_Index u, v;
    uint64_t weight;
} Edge;

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <graph_file>\n", argv[0]);
        return 1;
    }

    GrB_Info info;
    LAGraph_Init(NULL);

    FILE *file = fopen(argv[1], "r");
    if (!file) {
        perror("Error opening file");
        return 1;
    }

    char line[1024];
    GrB_Index num_nodes = 0;
    Edge *edges = NULL;
    size_t edge_capacity = 1024;
    size_t edge_count = 0;
    GrB_Index max_node = 0;

    edges = (Edge *)malloc(edge_capacity * sizeof(Edge));
    if (!edges) {
        perror("Memory allocation failed");
        fclose(file);
        return 1;
    }

    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '#' || line[0] == '\n')
            continue;

        GrB_Index source, neighbor, degree;
        uint64_t weight;
        char *ptr = line;

        if (sscanf(ptr, "%lu", &source) != 1) {
            fprintf(stderr, "Invalid line format: %s", line);
            fclose(file);
            free(edges);
            return 1;
        }

        if (source > max_node)
            max_node = source;

        while (*ptr && !isspace(*ptr))
            ptr++;
        while (*ptr && isspace(*ptr))
            ptr++;

        if (sscanf(ptr, "%lu", &degree) != 1) {
            fprintf(stderr, "Invalid degree in line: %s", line);
            fclose(file);
            free(edges);
            return 1;
        }

        while (*ptr && !isspace(*ptr))
            ptr++;
        while (*ptr && isspace(*ptr))
            ptr++;

        for (GrB_Index i = 0; i < degree; i++) {
            if (sscanf(ptr, "%lu %lu", &neighbor, &weight) != 2) {
                fprintf(stderr, "Invalid neighbor/weight in line: %s", line);
                fclose(file);
                free(edges);
                return 1;
            }

            if (neighbor > max_node)
                max_node = neighbor;

            if (source < neighbor) {
                if (edge_count >= edge_capacity) {
                    edge_capacity *= 2;
                    edges = (Edge *)realloc(edges, edge_capacity * sizeof(Edge));
                    if (!edges) {
                        perror("Memory reallocation failed");
                        fclose(file);
                        return 1;
                    }
                }

                edges[edge_count].u = source - 1;
                edges[edge_count].v = neighbor - 1;
                edges[edge_count].weight = weight;
                edge_count++;
            }

            while (*ptr && !isspace(*ptr))
                ptr++;
            while (*ptr && isspace(*ptr))
                ptr++;
            while (*ptr && !isspace(*ptr))
                ptr++;
            while (*ptr && isspace(*ptr))
                ptr++;
        }
    }
    fclose(file);

    num_nodes = max_node;

    if (num_nodes == 0 || edge_count == 0) {
        fprintf(stderr, "No graph data found in file\n");
        free(edges);
        return 1;
    }

    GrB_Matrix A;
    GrB_Matrix_new(&A, GrB_UINT64, num_nodes, num_nodes);

    for (size_t i = 0; i < edge_count; i++) {
        GrB_Matrix_setElement_UINT64(A, edges[i].weight, edges[i].u, edges[i].v);
        GrB_Matrix_setElement_UINT64(A, edges[i].weight, edges[i].v, edges[i].u);
    }

    GrB_Matrix msf;
    char msg[LAGRAPH_MSG_LEN];
    int status = LAGraph_msf(&msf, A, false, msg);

    if (status != GrB_SUCCESS) {
        printf("Error computing MSF: %s\n", msg);
        free(edges);
        GrB_free(&A);
        return 1;
    }

    GrB_Index nvals;
    GrB_Matrix_nvals(&nvals, msf);

    GrB_Index *I = (GrB_Index *)malloc(nvals * sizeof(GrB_Index));
    GrB_Index *J = (GrB_Index *)malloc(nvals * sizeof(GrB_Index));
    uint64_t *X = (uint64_t *)malloc(nvals * sizeof(uint64_t));

    GrB_Matrix_extractTuples_UINT64(I, J, X, &nvals, msf);

    for (GrB_Index k = 0; k < nvals; k++) {
        printf("%lu %lu %lu\n", I[k] + 1, J[k] + 1, X[k]);
    }

    free(I);
    free(J);
    free(X);
    free(edges);
    GrB_free(&A);
    GrB_free(&msf);

    return 0;
}