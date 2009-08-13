
  srand(10);
  for (int i=0 ; i<NUM ; i++) {
    double radius=1.0*rand()/RAND_MAX*5;
    double theta=1.0*rand()/RAND_MAX*2*M_PI; 
    r[i][0] = radius*cos(theta)+10;
    r[i][1] = radius*sin(theta)+10;
    r[i][2] = 0;
    /*
    r[i][0] = radius*1.0e-2+10;
    r[i][1] = radius*10+20;
    r[i][2] = 0;
    */

    radius=1.0*rand()/RAND_MAX*30;
    //theta=1.0*rand()/RAND_MAX*2*M_PI; 
    v[i][0] = radius*cos(theta+0.5*M_PI);
    v[i][1] = radius*sin(theta+0.5*M_PI);
    v[i][2] = 0.0;
/*
    v[i][0] = (1.0*rand()/RAND_MAX-0.5);
    v[i][1] = 0.0;
    v[i][2] = 0.0;
*/
/*
    for (int j=0 ; j<3 ; j++) {
      v[i][j] = 0.0;
      //v[i][j] = rand()/1e10;
    }
*/
    m[i] = 1.0;
  }


