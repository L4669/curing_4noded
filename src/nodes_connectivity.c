#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fem2d.h"

void node_connectivity(struct node_list *node_arr, struct ele_list *ele_arr, char *fname)
{
    char *node_line = "*node(";
    char *ele_line = "*quad4(";
    char *token;
    char temp[512];
    int node_inc=0;
    int ele_inc=0;

    FILE *fp = fopen(fname,"r");
    if(fp == NULL)
    {
        perror("In nodes_connectivity.c: File open failed\n");
        exit(0);
    }

    while(fgets(temp, 512, fp) != NULL)
    {
        if((strstr(temp, node_line)) != NULL)
        {
            token = strtok(temp, "(");
            token = strtok(NULL, "(");
            token = strtok(token, ",");
            node_arr[node_inc].node=atoi(token);
            token = strtok(NULL, ",");
            node_arr[node_inc].x=atof(token);
            token = strtok(NULL, ",");
            node_arr[node_inc].y=atof(token);
            token = strtok(NULL, ",");
            node_arr[node_inc].z=atof(token);
            node_inc++;
        }
        else if((strstr(temp, ele_line)) != NULL)
        {
            token = strtok(temp, "(");
            token = strtok(NULL, "(");
            token = strtok(token, ",");
            ele_arr[ele_inc].ele_no=atoi(token);
            token = strtok(NULL, ",");
            token = strtok(NULL, ",");
            ele_arr[ele_inc].ele_node_no[0]=atoi(token);
            token = strtok(NULL, ",");
            ele_arr[ele_inc].ele_node_no[1]=atoi(token);
            token = strtok(NULL, ",");
            ele_arr[ele_inc].ele_node_no[2]=atoi(token);
            token = strtok(NULL, ",");
            ele_arr[ele_inc].ele_node_no[3]=atoi(token);
            ele_inc++;
        }
    }

    fclose(fp);

}

void dirichlet_boundary_nodes(struct dirichlet_boundary_nodelist *db_array,
		char *fname, char *tag)
{
	FILE *fp;
	char buffer[512];
	char *token;
	int inc = 0, status = 0;

	fp = fopen(fname, "r");

    while(fgets(buffer, 512, fp) != NULL)
    {
        if((strstr(buffer, tag)) != NULL)
		{
			status = 1;
			continue;
		}
		else if((status == 1) && (strstr(buffer, "temperature") == NULL))
		{
			status = 0;
			break;
		}

		if(status == 1)
		{
			token = strtok(buffer, "(");
			token = strtok(NULL, ",");
			token = strtok(NULL, ",");
			token = strtok(NULL, ",");
			db_array[inc].node = atoi(token);
			token = strtok(NULL, ",");
			token = strtok(NULL, ",");
			db_array[inc].value = atof(token);
			inc++;
		}
	}

	fclose(fp);
}

void neumann_boundary_edges(struct neumann_boundary_edgelist *nb_array,
		struct ele_list *ele_arr, struct node_list *node_arr,
		char *fname, char *tag, int neu_edges)
{
	FILE *fp;
	char buffer[512];
	char *token;
	int line = 0, ele = 0, ele_no, node_local, edge_val;
	double *line_vector, dx, dy, dz, *product, vector_mag, 
		   dot_product, line_vector_mag;

	line_vector = (double *)malloc(sizeof(double)*3);
	product = (double *)malloc(sizeof(double)*3);

	fp = fopen(fname, "r");

	while(fgets(buffer, 512, fp) != NULL)
	{
		if(strstr(buffer, tag) != NULL)
		{
			break;
		}
	}

	for(line = 0; line < neu_edges; line++)
	{
		fgets(buffer, 512, fp);
		if(strstr(buffer, "pressure") != NULL)
		{
			token = strtok(buffer, "(");
			token = strtok(NULL, ",");
			token = strtok(NULL, ",");
			token = strtok(NULL, ",");
			nb_array[line].element = atoi(token);
			token = strtok(NULL, ",");
			token = strtok(NULL, ",");
			token = strtok(NULL, ",");
			nb_array[line].value = atof(token);
			token = strtok(NULL, ",");
			token = strtok(NULL, ",");
			nb_array[line].dcosine[0] = atof(token); 
			token = strtok(NULL, ",");
			nb_array[line].dcosine[1] = atof(token); 
			token = strtok(NULL, ",");
			nb_array[line].dcosine[2] = atof(token);
		}
		else
		{
			break;
		}
	}

	for(ele = 0; ele < neu_edges; ele++)
	{
		ele_no = nb_array[ele].element;
		nb_array[ele].edge_local = 0;
		for(node_local = 0; node_local < 4; node_local++)
		{
			dx = node_arr[ele_arr[ele_no - 1].ele_node_no[node_local] - 1].x
				- node_arr[ele_arr[ele_no - 1].ele_node_no[(node_local+1)%4] - 1].x;
			dy = node_arr[ele_arr[ele_no - 1].ele_node_no[node_local] - 1].y
				- node_arr[ele_arr[ele_no - 1].ele_node_no[(node_local+1)%4] - 1].y;
			dz = node_arr[ele_arr[ele_no - 1].ele_node_no[node_local] - 1].z
				- node_arr[ele_arr[ele_no - 1].ele_node_no[(node_local+1)%4] - 1].z;
			line_vector_mag = sqrt(dx*dx + dy*dy + dz*dz);
			//line_vector_mag = 1;
			line_vector[0] = -dx/line_vector_mag; 
			line_vector[1] = -dy/line_vector_mag; 
			line_vector[2] = -dz/line_vector_mag;

			dot_product = vector_dot_product(line_vector, nb_array[ele].dcosine);
			//printf("%d %f %f %f %f\n", ele+1, line_vector[0], line_vector[1], 
			//		line_vector[2], dot_product);

			if(fabs(dot_product) < 10e-2)
			{
				vector_cross_product(line_vector, nb_array[ele].dcosine, product);
				//printf("%d %f %f %f\n", ele+1, product[0], product[1], product[2]);
				if(product[2] > 0)
				{
					if(node_local+2 > 4)
					{
						edge_val = 10*(node_local+1) + ((node_local+2)%4);
					}
					else
					{
						edge_val = 10*(node_local+1) + (node_local+2);
					}
					nb_array[ele].edge_local = edge_val;
					nb_array[ele].length = sqrt(dx*dx + dy*dy + dz*dz);
				}
			}
		}
	}

	free(line_vector); free(product);
	//fclose to avoid sudden termination of program
	//causing termination of program when called inside the time loop
	fclose(fp);

}






