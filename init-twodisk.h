
  srand(10);
  for (int i=0 ; i<NUM/2 ; i++) {
    double radius=1.0*rand()/RAND_MAX*4;
    double theta=1.0*rand()/RAND_MAX*2*M_PI; 
    r[i][0] = radius*cos(theta)+14;
    r[i][1] = radius*sin(theta)+10;
    r[i][2] = 0;
    m[i] = 1.0;
  }

  for (int i=NUM/2 ; i<NUM ; i++) {
    double radius=1.0*rand()/RAND_MAX*4;
    double theta=1.0*rand()/RAND_MAX*2*M_PI; 
    r[i][0] = radius*cos(theta)+6;
    r[i][1] = radius*sin(theta)+10;
    r[i][2] = 0;
    m[i] = 1.0;
  }



  for (int i=0 ; i<NUM/2 ; i++) {
    double radius;
    double soft=5e-3;
    double mass=0.0;
    double length2=(r[i][0]-14)*(r[i][0]-14)+(r[i][1]-10)*(r[i][1]-10);
    double f[3];

    for (int k=0 ; k<3 ; k++)
      f[k] = 0.0;

    for (int j=0 ; j<NUM/2 ; j++) {
      if ( ((r[j][0]-14)*(r[j][0]-14)+(r[j][1]-10)*(r[j][1]-10)) < length2 ) {
        //mass += m[j];
        double r2;
        for (int k=0 ; k<3 ; k++)
          r2 += (r[j][k]-r[i][k])*(r[j][k]-r[i][k]);
        for (int k=0 ; k<3 ; k++)
          f[k] += (r[j][k]-r[i][k])/(r2+soft)/sqrt(r2+soft);
      }
    }
    //radius=sqrt(mass)/sqrt(length2+5e-2);
    radius=sqrt(sqrt(length2))*sqrt(sqrt(f[0]*f[0]+f[1]*f[1]+f[2]*f[2]));
    printf("debug: %lf\n",radius);
    v[i][0] = radius*(-((r[i][1]-10)/sqrt(length2)));
    v[i][1] = radius*((r[i][0]-14)/sqrt(length2))+6;
    v[i][2] = 0.0;
  }



  for (int i=NUM/2 ; i<NUM ; i++) {
    double radius;
    double soft=5e-3;
    double mass=0.0;
    double length2=(r[i][0]-6)*(r[i][0]-6)+(r[i][1]-10)*(r[i][1]-10);
    double f[3];

    for (int k=0 ; k<3 ; k++)
      f[k] = 0.0;

    for (int j=NUM/2 ; j<NUM ; j++) {
      if ( ((r[j][0]-6)*(r[j][0]-6)+(r[j][1]-10)*(r[j][1]-10)) < length2 ) {
        //mass += m[j];
        double r2;
        for (int k=0 ; k<3 ; k++)
          r2 += (r[j][k]-r[i][k])*(r[j][k]-r[i][k]);
        for (int k=0 ; k<3 ; k++)
          f[k] += (r[j][k]-r[i][k])/(r2+soft)/sqrt(r2+soft);

      }
    }

    //radius=sqrt(mass)/sqrt(length2+5e-2);
    radius=sqrt(sqrt(length2))*sqrt(sqrt(f[0]*f[0]+f[1]*f[1]+f[2]*f[2]));
    v[i][0] = radius*(-((r[i][1]-10)/sqrt(length2)));
    v[i][1] = radius*((r[i][0]-6)/sqrt(length2))-6;
    v[i][2] = 0.0;
  }

