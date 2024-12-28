#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string>

#include<windows.h>
#include<GL/glut.h>

#include<vector>

#include "bitmap_image.hpp"
#include "rayTracing_1705051.hpp"

#define PI 2*acos(0.0)
using namespace std;


#define MOVE_CONST 5
#define ROTATE_CONST ((0.5/180) * PI) // 0.5 degree

int windowWidth = 500;
int windowHeight = 500;
int captureCount  = 11;
double viewAngle = 80;
int recursion_level, pixels;
vector<Object*> objects;
vector<Light> lights;


int drawgrid;
int drawaxes;

// camera variables
Vector3D eye, upVec, rightVec, lookVec;

/* defining drawing functions */
void drawAxes()
{
	if(!drawaxes)
		return;

        glColor3f(1.0, 1.0, 1.0);
        glBegin(GL_LINES);{
            glVertex3f( 500,0,0);
            glVertex3f(-500,0,0);

            glVertex3f(0,500,0);
            glVertex3f(0, -500,0);

            glVertex3f(0,0,500);
            glVertex3f(0, 0,-500);
        }glEnd();

}

void drawSS()
{
    for (Object* o : objects) {
        o->draw();
    }

    for (Light l : lights) {
        l.draw();
    }
}


void capture() {
	cout << "Capturing bitmap image " <<pixels<< endl;
    int imageWidth, imageHeight, nearest, i, j, k;
	Vector3D topLeft, curPixel, temp;
    double t, tMin, du, dv;

	imageHeight = pixels;
	imageWidth = pixels;
    // initialize bitmap image
    bitmap_image image(pixels,pixels);

    // set background color
    for (i = 0; i < imageWidth; i++) {
        for (j = 0; j < imageHeight; j++) {
            image.set_pixel(i, j, 0, 0, 0);
        }
    }

    double planeDistance = (windowHeight * 0.5) / tan((viewAngle * 0.5) * (PI / 180));
	//print details
	cout << "Image Width: " << imageWidth << endl;
	cout << "Window Height: " << windowHeight << endl;
	cout << "Plane Distance: " << planeDistance << endl;
	cout << "View Angle: " << viewAngle << endl;
	cout << "Recursion Level: " << recursion_level << endl;
	cout << "Pixels: " << pixels << endl;
	cout << "Objects: " << objects.size() << endl;
	cout << "Lights: " << lights.size() << endl;
	temp = eye;
	temp = temp.add(lookVec.multiply(planeDistance));
	temp = temp.subtract(rightVec.multiply(windowWidth / 2.0));
	topLeft = temp.add(upVec.multiply(windowHeight / 2.0));


	//print calculation details
	cout << "Eye: " << eye << endl;

	cout << "Look Vector: " << lookVec << endl;
	cout << "Plane Distance: " << planeDistance << endl;
	temp = (lookVec.multiply(planeDistance));
	cout << "Look Vector * Plane Distance: " << temp << endl;

	cout << "Right Vector: " << rightVec << endl;
	cout << "Image Width: " << imageWidth << endl;
	temp = (rightVec.multiply(imageWidth / 2.0));
	cout << "Right Vector * Image Width/2: " << temp << endl;

	cout << "Up Vector: " << upVec << endl;
	cout << "Image Height: " << imageHeight << endl;
	temp = (upVec.multiply(imageHeight / 2.0));
	cout << "Up Vector * Image Height/2: " << temp << endl;

	cout << "Top Left: " << topLeft << endl;




    du = (double) windowWidth / imageWidth;
    dv = (double) windowHeight / imageHeight;
    // choose middle of the grid cell
	temp = rightVec.multiply(du / 2.0);
	temp = temp.subtract(upVec.multiply(dv / 2.0));
	topLeft = topLeft.add(temp);

	//print calculation details
	cout << "du: " << du << endl;
	cout << "Right Vector: " << rightVec << endl;
	temp = (rightVec.multiply(du / 2.0));
	cout << "Right Vector * du/2: " << temp << endl;

	cout << "dv: " << dv << endl;
	cout << "Up Vector: " << upVec << endl;
	temp = (upVec.multiply(dv / 2.0));
	cout << "Up Vector * dv/2: " << temp << endl;

	cout << "Top Left: " << topLeft << endl;



    for (i = 0; i < imageWidth; i++) {
        for (j = 0; j < imageHeight; j++) {
            nearest = -1;
            tMin = INFINITY;
            // calculate curPixel
			temp = rightVec.multiply(du * i);
			temp = temp.subtract(upVec.multiply(dv * j));
			curPixel = topLeft.add(temp);
            // cast ray
            Ray ray(eye, (curPixel.subtract(eye)));
            Color color;

            int obj_size = objects.size();
            for (k = 0; k < obj_size; k++) {
				bool is_print = false;
				//if(j==0 && i<3){is_print = true;}
				if(is_print){
					//print i and j
					cout<<"==================================================="<<endl;
					cout << "i: " << i << endl;
					cout << "j: " << j << endl;
					cout<<"==================================================="<<endl;
					cout<<"Input Color: "<< color <<endl;
				}
                t = objects[k]->intersectWithIllumination(ray, color, 0, is_print);
				if(is_print) cout<<"Output Color: "<< color <<endl;

				if (t < tMin) {
                    if(t > 0){
						nearest = k;
						tMin = t;
					}
                }
            }

            if (nearest != -1) {
				double tNear = objects[nearest]->intersectWithIllumination(ray, color, 1);
                tMin = tNear;
            }
            color.clip();

            image.set_pixel(i, j, (color.getR() * 255), (color.getG() * 255), (color.getB()) * 255);
        }
    }
    string outPath = "C:\\Users\\USER\\Desktop\\RayTracing\\output_";
    outPath += to_string(captureCount)+".bmp";
    captureCount++;
    image.save_image(outPath);
    image.clear();
	cout << "Finished Capturing bitmap image. Path: " <<outPath << endl;
}



/* defining listener functions */

void keyboardListener(unsigned char key, int x,int y){
    switch(key){
        case '0':
            capture();
            break;
        case '1':
            rightVec = rotateVector(rightVec, upVec, ROTATE_CONST);
            lookVec = rotateVector(lookVec, upVec, ROTATE_CONST);
            break;
        case '2':
            rightVec = rotateVector(rightVec, upVec, -ROTATE_CONST);
            lookVec = rotateVector(lookVec, upVec, -ROTATE_CONST);
            break;
        case '3':
            upVec = rotateVector(upVec, rightVec, ROTATE_CONST);
            lookVec = rotateVector(lookVec, rightVec, ROTATE_CONST);
            break;
        case '4':
            upVec = rotateVector(upVec, rightVec, -ROTATE_CONST);
            lookVec = rotateVector(lookVec, rightVec, -ROTATE_CONST);
            break;
        case '5':
            upVec = rotateVector(upVec, lookVec, -ROTATE_CONST);
            rightVec = rotateVector(rightVec, lookVec, -ROTATE_CONST);
            break;
        case '6':
            upVec = rotateVector(upVec, lookVec, ROTATE_CONST);
            rightVec = rotateVector(rightVec, lookVec, ROTATE_CONST);
            break;
        default:
            break;
    }
}

void specialKeyListener(int key, int x,int y){
	Vector3D t;
    switch(key){
        case GLUT_KEY_UP:		// upVec arrow key
			t = lookVec.multiply(MOVE_CONST);
            eye = eye.add(t);
            break;
        case GLUT_KEY_DOWN:		//down arrow key
            t = lookVec.multiply(MOVE_CONST);
            eye = eye.subtract(t);
            break;
        case GLUT_KEY_RIGHT:
            t = rightVec.multiply(MOVE_CONST);
            eye = eye.add(t);
            break;
        case GLUT_KEY_LEFT:
            t = rightVec.multiply(MOVE_CONST);
            eye = eye.subtract(t);
            break;
        case GLUT_KEY_PAGE_UP:
            t = upVec.multiply(MOVE_CONST);
            eye = eye.add(t);
            break;
        case GLUT_KEY_PAGE_DOWN:
            t = rightVec.multiply(MOVE_CONST);
            eye = eye.subtract(t);
            break;
        default:
            break;
    }
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
    switch(button){
        case GLUT_LEFT_BUTTON:
            if(state == GLUT_DOWN) {

            }
            break;
        case GLUT_RIGHT_BUTTON:
            if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
                drawaxes=1-drawaxes;
            }
            break;
        default:
            break;
    }
}


/* defining auxiliary functions */

void display(){

    //clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0,0,0,0);	//color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    /********************
    / set-upVec camera here
    ********************/
    //load the correct matrix -- MODEL-VIEW matrix
    glMatrixMode(GL_MODELVIEW);

    //initialize the matrix
    glLoadIdentity();

    //now give three info
    //1. where is the camera (viewer)?
    //2. where is the camera looking?
    //3. Which direction is the camera's UP direction?

    gluLookAt(eye.getX(), eye.getY(), eye.getZ(), eye.getX() + lookVec.getX(), eye.getY() + lookVec.getY(), eye.getZ() + lookVec.getZ(), upVec.getX(), upVec.getY(), upVec.getZ());


    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    /****************************
    / Add your objects from here
    ****************************/
    //add objects

    drawAxes();

    drawSS();

    //ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}

void animate() {
	/* codes for any changes in Models, Camera */
	glutPostRedisplay();
}

void init() {
    //codes for initialization
    drawgrid=0;
    drawaxes=1;

	double val = 1.0/(sqrt(2.0));

    // initialize camera variables
	eye.setVector(100, 100, 0);

	upVec.setVector(0, 0, 1);
    upVec.normalize();

	rightVec.setVector(-val, val, 0);
    rightVec.normalize();

	lookVec.setVector(-val, -val, 0);
    lookVec.normalize();

    //clear the screen
    glClearColor(0,0,0,0);

    /************************
    / set-upVec projection here
    ************************/
    //load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    //initialize the matrix
    glLoadIdentity();

    //give PERSPECTIVE parameters
    gluPerspective(viewAngle,	1,	1,	1000.0);
    //field of view in the Y (vertically)
    //aspect ratio that determines the field of view in the X direction (horizontally)
    //near distance
    //far distance
}

void clearObjects() {
    /* clearing objects vector */
    objects.clear();
}

void clearLights() {
    /* clearing lights vector */
    lights.clear();
}

void loadData() {
    ifstream input;

    /* preparing input for extracting values from input file */

    input.open("C:\\Users\\USER\\Desktop\\RayTracing\\scene.txt");

	if (!input) {
        cout << "Unable to open file scene.txt" << endl;
        exit(1);
    }

    /* extracting recursion level,image pixel dimension, & object count from input file */
    int object_count;
	input >> recursion_level;
	input >> pixels;
	input  >> object_count;



    string objectShape;

    Object* object = NULL;

    for(int i=0; i<object_count; i++) {
        input >> objectShape;

        if(objectShape == "sphere") {
            Vector3D center;
            double radius;

            input >> center;
            input >> radius;

            object = new Sphere(center, radius);
        } else if(objectShape == "triangle") {
            Vector3D a, b, c;

            input >> a >> b >> c;

            object = new Triangle(a, b, c);
        } else if(objectShape == "general") {
            Vector3D cubeReferencePoint;
            double length, width, height;

			double A, B, C, D, E, F, G, H, I, J;
            input >> A >> B >> C >> D >> E >> F >> G >> H >> I >> J;

            input >> cubeReferencePoint;
            input >> length;
			input >> width;
			input >> height;

			object = new GeneralQuadricSurface(A, B, C, D, E, F, G, H, I, J, length, width, height, cubeReferencePoint);
        }
		else {
            cout << objectShape;
			cout << ": invalid object shape found" << endl;
            break;
        }

        Color color;
		input >> color;
		object->setColor(color);

        ReflectionCoefficients reflectionCoefficient;
		input >> reflectionCoefficient;
		object->setCoefficients(reflectionCoefficient);

        int shininess;
		input >> shininess;
		object->setShine(shininess);

        objects.push_back(object);
    }
    object = NULL;

    /* extracting point light sources information from input file */
    int lightsCount;
	input >> lightsCount;
	int i;

    for(i=0; i<lightsCount; i++) {
        Vector3D position;
		input >> position;

		Color color;
        input >> color;

        lights.push_back(Light(position, color));
    }

	/* extracting spot light sources information from input file */
	input >> lightsCount;

    for(i=0; i<lightsCount; i++) {
        Vector3D position;
		input >> position;

		Color color;
        input >> color;

		Vector3D direction;
		input >> direction;

		double cutoffAngle;
		input >> cutoffAngle;

		Light l(position, color);
		l.setSpotLight(direction, cutoffAngle);

        lights.push_back(l);
    }
    input.close();

    /* creating a floor object and pushing it to objects vector */
    object = new Floor(1000.0, 20.0);
    object->setCoefficients(0.3, 0.3, 0.3, 0.3);
    object->setShine(30);

    objects.push_back(object);

	//printing objects information
	int obj_size = objects.size();
	for(int i=0; i<obj_size; i++) {
		cout << "object " << i << ": " << objects[i]->getLength() << endl;
	}
    object = NULL;
}


/* defining main function */
int main(int argc, char **argv){
    loadData();

    glutInit(&argc,argv);
    glutInitWindowSize(windowWidth, windowHeight);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

    glutCreateWindow("My OpenGL Program");

    init();

    glEnable(GL_DEPTH_TEST);	//enable Depth Testing

    glutDisplayFunc(display);	//display callback function
    glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occurring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    glutMainLoop();		//The main loop of OpenGL

    clearObjects();
	clearLights();
    return 0;
}
