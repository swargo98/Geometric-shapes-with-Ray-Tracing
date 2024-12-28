#ifndef RAYTRACING_1705051_HPP_INCLUDED
#define RAYTRACING_1705051_HPP_INCLUDED

#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include<limits>
#include <cstdio>
#include <GL/glut.h>


using namespace std;

#define PI 2*acos(0.0)
#define INF numeric_limits<double>::infinity()
#define EPSILON 0.0000001

class Object;
class Light;


extern vector<Object*> objects;
extern vector<Light> lights;
extern int recursion_level;

class Color {
    double clipValue(double value){
        if(value > 1.0)
            value = 1.0;
        if(value < 0.0)
            value = 0.0;
        return value;
    }
    double r, g, b;
public:
    Color() {
        r = 0;
        g = 0;
        b = 0;
        }
    Color(double r, double g, double b){
        this->r = r;
        this->g = g;
        this->b = b;
    }
    double getR(){
        return r;
    }
    double getG(){
        return g;
    }
    double getB(){
        return b;
    }
    void setR(double r){
        this->r = r;
    }
    void setG(double g){
        this->g = g;
    }
    void setB(double b){
        this->b = b;
    }
    void setColor(double r, double g, double b){
        this->r = r;
        this->g = g;
        this->b = b;
    }
    Color add(Color c){
        return Color(r+c.r, g+c.g, b+c.b);
    }
    Color multiply(double k){
        return Color(r*k, g*k, b*k);
    }
    Color multiply(Color c){
        return Color(r*c.r, g*c.g, b*c.b);
    }
    void clip(){
        r = clipValue(r);
        g = clipValue(g);
        b = clipValue(b);
    }

    friend ifstream& operator>>(ifstream&, Color&);
    friend ostream& operator<<(ostream&, Color&);
};

ifstream &operator>>(ifstream &in, Color &c){
    in >> c.r >> c.g >> c.b;
    return in;
}

ostream &operator<<(ostream &out, Color &c){
    out << '[' << c.r << ", " << c.g << ", " << c.b  << ']' << endl;
    return out;
}

class Vector3D {
    double x;
    double y;
    double z;

public:
    Vector3D() {
        x = 0.0;
        y = 0.0;
        z = 0.0;
    }

    Vector3D(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    double getX(){
        return x;
    }

    double getY(){
        return y;
    }

    double getZ(){
        return z;
    }

    void setX(double x){
        this->x = x;
    }

    void setY(double y){
        this->y = y;
    }

    void setZ(double z){
        this->z = z;
    }

    void setVector(double x, double y, double z){
        this->x = x;
        this->y = y;
        this->z = z;
    }

    void normalize(){
        double length;
        length = sqrt(x*x + y*y + z*z);
        x = x / length;
        y = y / length;
        z = z / length;
    }
    double computeDistanceBetween(Vector3D v){
        return sqrt((x-v.getX())*(x-v.getX()) + (y-v.getY())*(y-v.getY()) + (z-v.getZ())*(z-v.getZ()));
    }

    Vector3D add(Vector3D v){
        return Vector3D(x+v.getX(), y+v.getY(), z+v.getZ());
    }
    Vector3D subtract(Vector3D v){
        return Vector3D(x-v.getX(), y-v.getY(), z-v.getZ());
    }
    Vector3D multiply(double k){
        return Vector3D(x*k, y*k, z*k);
    }
    Vector3D crossProduct(Vector3D v){
        return Vector3D(y*v.getZ() - z*v.getY(), z*v.getX() - x*v.getZ(), x*v.getY() - y*v.getX());
    }
    double dotProduct(Vector3D v){
        return x*v.getX() + y*v.getY() + z*v.getZ();
    }

    friend ifstream& operator>>(ifstream&, Vector3D&);
    friend ostream& operator<<(ostream&, Vector3D&);
};

ifstream &operator>>(ifstream &in, Vector3D &v){
    in >> v.x >> v.y >> v.z;
    return in;
}

ostream &operator<<(ostream &out, Vector3D &v){
    out << '[' << v.x << ", " << v.y << ", " << v.z  << ']' << endl;
    return out;
}

Vector3D rotateVector(Vector3D a, Vector3D b, double angle){
    Vector3D v = (a.multiply(cos(angle)));
    Vector3D h = (a.crossProduct(b));
    h = h.multiply(sin(angle));
    Vector3D result = v.add(h);
    result.normalize();
    return result;
}

//compute the angle between two vectors
double computeAngle(Vector3D a, Vector3D b){
    double dot = a.dotProduct(b);
    double a_length = a.computeDistanceBetween(Vector3D(0, 0, 0));
    double b_length = b.computeDistanceBetween(Vector3D(0, 0, 0));
    double angle = acos(dot/(a_length*b_length));
    angle = angle * 180 / PI;
    return angle;
}
class Ray {
    Vector3D origin;
    Vector3D direction;
public:
    Ray() {
        origin = Vector3D(0.0, 0.0, 0.0);
        direction = Vector3D(0.0, 0.0, 0.0);
    }

    Ray(Vector3D origin, Vector3D direction) {
        this->origin = origin;
        this->direction = direction;
        this->direction.normalize();
    }

    Vector3D getOrigin() {
        return origin;
    }

    Vector3D getDirection() {
        return direction;
    }
};

class Light {
    Vector3D light_pos;
    Color color;
    bool is_SpotLight;
    Vector3D spotDirection;
    double spotCutoff;

    double radius;
    int segments;
    int stacks;
public:
    Light(){
        light_pos = Vector3D(0.0, 0.0, 0.0);
        color = Color(0.0, 0.0, 0.0);
        radius = 0.0;
        segments = 0;
        stacks = 0;
        is_SpotLight = false;
        spotDirection = Vector3D(0.0, 0.0, 0.0);
        spotCutoff = 360.0;
    }


    Light(Vector3D l_pos, Color c, double r = 1.0, int seg = 12, int stck = 4){
        light_pos = l_pos;
        color = c;
        radius = r;
        segments = seg;
        stacks = stck;
        is_SpotLight = false;
        spotDirection = Vector3D(0.0, 0.0, 0.0);
        spotCutoff = 360.0;
    }

    void setSpotLight(Vector3D dir, double cutoff){
        is_SpotLight = true;
        spotDirection = dir;
        spotCutoff = cutoff;
        cout<<"Spotlight is set"<<endl;
        cout<<"Spotlight direction: "<<spotDirection<<endl;
        cout<<"Spotlight cutoff: "<<spotCutoff<<endl;
    }

    bool isSpotLight(){
        return is_SpotLight;
    }

    Vector3D getSpotDirection(){
        return spotDirection;
    }

    double getSpotCutoff(){
        return spotCutoff;
    }

    Color getColor(){
        return color;
    }

    Vector3D getLightPos(){
        return light_pos;
    }



    void draw(){
        stacks = stacks + 1;
        segments = segments + 1;
        Vector3D points[stacks][segments];
        int i,j;
        double height, _radius;

        /* generating points: segments = segments in plane; stacks = segments in hemisphere */
        for(i=0; i<stacks; i++) {
            height = radius*sin(((double)i/(double)(stacks-1))*(PI/2));
            _radius = radius*cos(((double)i/(double)(stacks-1))*(PI/2));

            for(j=0; j<segments; j++) {
                points[i][j] = Vector3D(_radius*cos(((double)j/(double)(segments-1))*2*PI), _radius*sin(((double)j/(double)(segments-1))*2*PI), height);
            }
        }

        /* drawing quads using generated points */
        glColor3f(color.getR(), color.getG(), color.getB());

        for(i=0; i<(stacks-1); i++) {
            for(j=0; j<(segments-1); j++) {
                glBegin(GL_QUADS);
                {
                    /* upper hemisphere */
                    glVertex3f(
                        (light_pos.add(points[i][j])).getX(),
                        (light_pos.add(points[i][j])).getY(),
                        (light_pos.add(points[i][j])).getZ()
                        );

                    glVertex3f(
                        (light_pos.add(points[i][j+1])).getX(),
                        (light_pos.add(points[i][j+1])).getY(),
                        (light_pos.add(points[i][j+1])).getZ()
                        );
                    glVertex3f(
                        (light_pos.add(points[i+1][j+1])).getX(),
                        (light_pos.add(points[i+1][j+1])).getY(),
                        (light_pos.add(points[i+1][j+1])).getZ()
                        );
                    glVertex3f(
                        (light_pos.add(points[i+1][j])).getX(),
                        (light_pos.add(points[i+1][j])).getY(),
                        (light_pos.add(points[i+1][j])).getZ()
                        );

                    /* lower hemisphere */
                    glVertex3f(
                        (light_pos.add(points[i][j])).getX(),
                        (light_pos.add(points[i][j])).getY(),
                        (light_pos.subtract(points[i][j])).getZ()
                        );
                    glVertex3f(
                        (light_pos.add(points[i][j+1])).getX(),
                        (light_pos.add(points[i][j+1])).getY(),
                        (light_pos.subtract(points[i][j+1])).getZ()
                        );
                    glVertex3f(
                        (light_pos.add(points[i+1][j+1])).getX(),
                        (light_pos.add(points[i+1][j+1])).getY(),
                        (light_pos.subtract(points[i+1][j+1])).getZ()
                        );
                    glVertex3f(
                        (light_pos.add(points[i+1][j])).getX(),
                        (light_pos.add(points[i+1][j])).getY(),
                        (light_pos.subtract(points[i+1][j])).getZ()
                        );
                }
                glEnd();
            }
        }
    }
};

class ReflectionCoefficients{
    double ambientReflectionCoefficient;
    double diffuseReflectionCoefficient;
    double specularReflectionCoefficient;
    double recursiveReflectionCoefficient;
    public:
    ReflectionCoefficients(){
        ambientReflectionCoefficient = 0.0;
        diffuseReflectionCoefficient = 0.0;
        specularReflectionCoefficient = 0.0;
        recursiveReflectionCoefficient = 0.0;
    }

    ReflectionCoefficients(double ambientReflectionCoefficient, double diffuseReflectionCoefficient, double specularReflectionCoefficient, double recursiveReflectionCoefficient){
        this->ambientReflectionCoefficient = ambientReflectionCoefficient;
        this->diffuseReflectionCoefficient = diffuseReflectionCoefficient;
        this->specularReflectionCoefficient = specularReflectionCoefficient;
        this->recursiveReflectionCoefficient = recursiveReflectionCoefficient;
    }

    double getAmbientReflectionCoefficient(){
        return ambientReflectionCoefficient;
    }
    double getDiffuseReflectionCoefficient(){
        return diffuseReflectionCoefficient;
    }
    double getSpecularReflectionCoefficient(){
        return specularReflectionCoefficient;
    }
    double getRecursiveReflectionCoefficient(){
        return recursiveReflectionCoefficient;
    }
    void setAmbientReflectionCoefficient(double ambientReflectionCoefficient){
        this->ambientReflectionCoefficient = ambientReflectionCoefficient;
    }
    void setDiffuseReflectionCoefficient(double diffuseReflectionCoefficient){
        this->diffuseReflectionCoefficient = diffuseReflectionCoefficient;
    }
    void setSpecularReflectionCoefficient(double specularReflectionCoefficient){
        this->specularReflectionCoefficient = specularReflectionCoefficient;
    }
    void setRecursiveReflectionCoefficient(double recursiveReflectionCoefficient){
        this->recursiveReflectionCoefficient = recursiveReflectionCoefficient;
    }

    friend ifstream& operator>>(ifstream&, ReflectionCoefficients&);
    friend ostream& operator<<(ostream&, ReflectionCoefficients&);
};

ifstream& operator>>(ifstream& in, ReflectionCoefficients& ref){
    in >> ref.ambientReflectionCoefficient;
    in >> ref.diffuseReflectionCoefficient;
    in >> ref.specularReflectionCoefficient;
    in >> ref.recursiveReflectionCoefficient;
    return in;
}

ostream& operator<<(ostream& out, ReflectionCoefficients& ref){
    out << '[';
    out << ref.ambientReflectionCoefficient << ", ";
    out << ref.diffuseReflectionCoefficient << ", ";
    out << ref.specularReflectionCoefficient << ", ";
    out << ref.recursiveReflectionCoefficient << ", ";
    out << ']';
    return out;
}

class Object {
protected:
    Vector3D reference_point;
    double height, width, length;
    Color color;
    ReflectionCoefficients coefficients; // reflection coefficients
    int shine; // exponent term of specular component

public:
    Object()= default;

    double getLength(){
        return length;
    }

    virtual void draw(){}

    virtual double intersect(Ray r, Color clr, int level, bool is_print = false){
        return -1.0;
    }

    virtual Vector3D getNormalAt(Vector3D intersectionPoint) {
        return Vector3D(0.0, 0.0, 0.0);
    }

    virtual Color getColorAt(Vector3D intersectionPoint) {
        return color;
    }

    double intersectWithIllumination(Ray& r, Color& clr, int level, bool is_print = false){
        double tmin;
        tmin = this->intersect(r, clr, level, is_print);
        //print details of the intersection
        if(is_print){
            cout << "tmin: " << tmin << endl;
            cout << "clr: " << clr << endl;
            //print ray details
            Vector3D temp = r.getOrigin();
            cout << "ray.getOrigin(): " << temp << endl;
            temp = r.getDirection();
            cout << "ray.getDirection(): " << temp << endl;
        }

        if(level == 0) return tmin;

        // illumination with phong lighting model
        Vector3D ro, rd, intersectionPoint, normal;
        ro = r.getOrigin(); // origin
        rd = r.getDirection(); // direction
        intersectionPoint = ro.add( rd.multiply(tmin) );

        clr = getColorAt(intersectionPoint).multiply(coefficients.getAmbientReflectionCoefficient()); // ambient
        clr.clip();

        normal = getNormalAt(intersectionPoint);
        normal.normalize();

        for (Light l : lights) {
            Vector3D lightDir, lightPos;
            Ray lightRay;

            // cast ray from light source
            lightDir = l.getLightPos();
            lightDir = lightDir.subtract(intersectionPoint);
            lightDir.normalize();

            // prevent self intersection by slightly moving start in the direction of light
            lightPos = intersectionPoint.add(lightDir.multiply(0.0000000001));
            lightRay = Ray(lightPos, lightDir);

            if(l.isSpotLight()){
                Vector3D t = intersectionPoint.subtract(l.getLightPos());
                double angle = computeAngle(t, l.getSpotDirection());
                if(angle > l.getSpotCutoff()) continue;
            }


            // check if any other object is present between this object & light source
            bool inShadow;
            Color temp;
            double t, tMinActual;

            inShadow = false;
            tMinActual = INFINITY;
            int numObjects = objects.size();
            for (int i = 0; i < numObjects; i++) {
                t = objects[i]->intersectWithIllumination(lightRay, temp, 0);
                if (t > 0) {
                    if(t < tMinActual) tMinActual = t;
                }
            }
            if (tmin <= tMinActual) {
                inShadow = inShadow;
            }
            else {
                inShadow = true;
            }

            // compute diffuse and specular components
            if (!inShadow) {
                double lambert, phong;
                lambert = max(normal.dotProduct(lightDir), 0.0);

                // R = 2(L.N)N – L
                Vector3D R;
                R = normal.multiply(2.0);
                R = R.multiply(normal.dotProduct(lightDir));
                R = R.subtract(lightDir);
                R.normalize();

                phong = max(pow((rd.dotProduct(R)), shine), 0.0);
                Color temp = l.getColor().multiply(coefficients.getDiffuseReflectionCoefficient() * lambert);
                temp = temp.multiply(getColorAt(intersectionPoint));
                clr = clr.add(temp);
                clr.clip();
                temp = l.getColor().multiply(coefficients.getSpecularReflectionCoefficient() * phong);
                clr = clr.add(temp);
                clr.clip();
            }

        }

        // recursive reflection
        if (level >= recursion_level) return tmin;
        // construct reflected ray from intersection point
        Vector3D rayReflectedDir, rayReflectedOrigin, temp_v;

        temp_v = normal.multiply(2.0);
        temp_v = temp_v.multiply(normal.dotProduct(rd));

        rayReflectedDir = rd.subtract(temp_v);
        rayReflectedDir.normalize();
        rayReflectedOrigin = intersectionPoint.add(rayReflectedDir.multiply(0.0000000001));
        Ray rayReflected(rayReflectedOrigin, rayReflectedDir);

        // find nearest intersecting object
        Color colorReflected;
        int nearest = -1;
        double tReflected, tReflectedMin;

        tReflectedMin = INFINITY;
        int objCount = objects.size();
        for (int k = 0; k < objCount; k++) {
            tReflected = objects[k]->intersectWithIllumination(rayReflected, colorReflected, 0);
            if (tReflected > 0) {
                if(tReflected < tReflectedMin){
                    nearest = k;
                    tReflectedMin = tReflected;
                } 
            }
        }

        if (nearest != -1) {
            tReflectedMin = objects[nearest]->intersectWithIllumination(rayReflected, colorReflected, level + 1);

            Color temp = colorReflected.multiply(coefficients.getRecursiveReflectionCoefficient());
            clr = clr.add(temp);
            clr.clip();
        }
        return tmin;
    }

    void setColor(Color c){
        color = c;
    }

    void setColor(double c1, double c2, double c3){
        color.setColor(c1, c2, c3);
    }
    void setShine(int s){
        shine = s;
    }
    void setCoefficients(double c1, double c2, double c3, double c4){
        coefficients.setAmbientReflectionCoefficient(c1);
        coefficients.setDiffuseReflectionCoefficient(c2);
        coefficients.setSpecularReflectionCoefficient(c3);
        coefficients.setRecursiveReflectionCoefficient(c4);
    }

    void setCoefficients(ReflectionCoefficients c){
        coefficients = c;
    }
};

class Sphere : public Object {
public:
    Sphere(Vector3D center, double radius) {
        reference_point = center;
        length = radius;
    }

    void draw() override{
        glTranslatef(reference_point.getX(), reference_point.getY(), reference_point.getZ());
        int stacks, slices;
        stacks = 20;
        slices = 24;
        Vector3D points[100][100];
        double radius;
        radius = length;
        //generate points
        for(int i=0;i<(stacks+1);i++)
        {
            double h,r, angle;
            angle = ((double)i/(double)stacks)*(PI/2);
            h=radius*sin(angle);
            r=radius*cos(angle);
            for(int j=0;j<(slices+1);j++)
            {
                angle = ((double)j/(double)slices)*2*PI;
                points[i][j].setX(r*cos(angle));
                points[i][j].setY(r*sin(angle));
                points[i][j].setZ(h);
            }
        }
        //draw quads using generated points
        glColor3f(color.getR(), color.getG(), color.getB());
        for(int i=0;i<stacks;i++)
        {
            for(int j=0;j<slices;j++)
            {
                glBegin(GL_QUADS);{
                    //upper hemisphere
                    glVertex3f(points[i][j].getX(),points[i][j].getY(),points[i][j].getZ());
                    glVertex3f(points[i][j+1].getX(),points[i][j+1].getY(),points[i][j+1].getZ());
                    glVertex3f(points[i+1][j+1].getX(),points[i+1][j+1].getY(),points[i+1][j+1].getZ());
                    glVertex3f(points[i+1][j].getX(),points[i+1][j].getY(),points[i+1][j].getZ());
                    //lower hemisphere
                    glVertex3f(points[i][j].getX(),points[i][j].getY(),-points[i][j].getZ());
                    glVertex3f(points[i][j+1].getX(),points[i][j+1].getY(),-points[i][j+1].getZ());
                    glVertex3f(points[i+1][j+1].getX(),points[i+1][j+1].getY(),-points[i+1][j+1].getZ());
                    glVertex3f(points[i+1][j].getX(),points[i+1][j].getY(),-points[i+1][j].getZ());
                }glEnd();
            }
        }
        glTranslatef(-reference_point.getX(), -reference_point.getY(), -reference_point.getZ());
    }
    double intersect (Ray r, Color clr, int level, bool is_print = false) override{
        //cout<<"Sphere intersect"<<endl;
        Vector3D ro = (r.getOrigin()).subtract(reference_point); // origin
        Vector3D rd = r.getDirection(); // direction
        double radius, a, b, c, d, temp, t1, t2;
        radius = length;
        a = 1;
        b = (rd.dotProduct(ro)) * 2;
        c = ((ro.dotProduct(ro))) - (radius * radius);
        temp = b * b - 4 * a * c;
        if(is_print){
            cout << "ro: " << ro.getX() << " " << ro.getY() << " " << ro.getZ() << endl;
            cout << "rd: " << rd.getX() << " " << rd.getY() << " " << rd.getZ() << endl;
            cout << "a: " << a << endl;
            cout << "b: " << b << endl;
            cout << "c: " << c << endl;
            cout << "temp: " << temp << endl;
        }
        if (temp < 0) return -1;
        d = sqrt(temp);
        t1 = (- b - d) / (2 * a);
        t2 = (- b + d) / (2 * a);

        if(is_print){
            cout << "t1: " << t1 << endl;
            cout << "t2: " << t2 << endl;
        }

        if (t1 < 0){
            if(t2 < 0) return -1;
        }
        else if (t1 > 0){
            return t1;
        }
        else if (t2 > 0){
            return t2;
        }
        return -1;
    }
    Vector3D getNormalAt (Vector3D intersectionPoint) override {
        Vector3D n = intersectionPoint.subtract(reference_point);
        n.normalize();
        return n;
    }
};


class Triangle : public Object {
    Vector3D v1, v2, v3;
public:
    Triangle(Vector3D p1, Vector3D p2, Vector3D p3) {
        v1 = p1;
        v2 = p2;
        v3 = p3;
    }

    void draw() override{
        glBegin(GL_TRIANGLES);{
            glColor3f(color.getR(), color.getG(), color.getB());
            glVertex3f(v1.getX(), v1.getY(), v1.getZ());
            glVertex3f(v2.getX(), v2.getY(), v2.getZ());
            glVertex3f(v3.getX(), v3.getY(), v3.getZ());
        }glEnd();
    }
    double intersect (Ray r, Color clr, int level, bool is_print = false) override{
        Vector3D ro, rd, edge1, edge2, h, s, q;
        double a, f, t, u, v;

        ro = r.getOrigin();
        rd = r.getDirection();

        edge1 = v2.subtract(v1);
        edge2 = v3.subtract(v1);

        h = rd.crossProduct(edge2);
        a = edge1.dotProduct(h);

        if (a > -EPSILON && a < EPSILON) {
            return -1; // ray is parallel to the triangle
        }

        f = 1.0 / a;
        s = ro.subtract(v1);
        u = f * s.dotProduct(h);

        if (u < 0.0 || u > 1.0) {
            return -1;
        }

        q = s.crossProduct(edge1);
        v = f * rd.dotProduct(q);

        if (v < 0.0 || (u+v) > 1.0) {
            return -1;
        }

        t = f * edge2.dotProduct(q);
        if (t > EPSILON) {
            return t; // ray intersection
        } else {
            return -1; // line intersection
        }
    }
    Vector3D getNormalAt (Vector3D intersectionPoint) override{
        Vector3D edge1, edge2, n;

        edge1 = v2.subtract(v1);
        edge2 = v3.subtract(v1);

        n = edge1.crossProduct(edge2);
        n.normalize();
        return n;
    }
};

class GeneralQuadricSurface : public Object {
    double A, B, C, D, E, F, G, H, I, J;

    bool withinReferenceCube(Vector3D p){
        bool within;
        within = true;
        
        if (height != 0) {
            if (p.getZ() < reference_point.getZ()) within = false;
            if (p.getZ() > reference_point.getZ() + height) within = false;
            
        }
        if (length != 0) {
            if (p.getX() < reference_point.getX()) within = false;
            if (p.getX() > reference_point.getX() + length) within = false;
        }
        if (width != 0) {
            if (p.getY() < reference_point.getY()) within = false;
            if (p.getY() > reference_point.getY() + width) within = false;
        }
        
        return within;
    }
public:
    GeneralQuadricSurface(double a, double b, double c, double d, double e, double f, double g, double h, double i, double j, double length, double width, double height, Vector3D ref) {
        A = a;
        B = b;
        C = c;
        D = d;
        E = e;
        F = f;
        G = g;
        H = h;
        I = i;
        J = j;
        this->length = length;
        this->width = width;
        this->height = height;
        reference_point = ref;
    }
    double intersect (Ray r, Color clr, int level, bool is_print = false) override{
        Vector3D ro, rd, p1, p2;
        double a, b, c, d, temp, t, t1, t2;

        ro = r.getOrigin();
        rd = r.getDirection();

        a = A * rd.getX() * rd.getX() + B * rd.getY() * rd.getY() + C * rd.getZ() * rd.getZ();
        a += D * rd.getX() * rd.getY() + E * rd.getX() * rd.getZ() + F * rd.getY() * rd.getZ();

        b = 2 * (A * rd.getX() * ro.getX() + B * rd.getY() * ro.getY() + C * rd.getZ() * ro.getZ());
        b += D * (rd.getX() * ro.getY() + rd.getY() * ro.getX()) + E * (rd.getX() * ro.getZ() + rd.getZ() * ro.getX());
        b += F * (rd.getY() * ro.getZ() + rd.getZ() * ro.getY()) + G * rd.getX() + H * rd.getY() + I * rd.getZ();

        c = A * ro.getX() * ro.getX() + B * ro.getY() * ro.getY() + C * ro.getZ() * ro.getZ();
        c += D * (ro.getX() * ro.getY()) + E * (ro.getX() * ro.getZ()) + F * (ro.getY() * ro.getZ());
        c += G * ro.getX() + H * ro.getY() + I * ro.getZ() + J;


        temp = b * b;
        temp -= 4 * a * c;
        if (temp < 0){
            return -1;
        }

        d = sqrt(temp);
        t1 = (- b - d) / (2 * a);
        t2 = (- b + d) / (2 * a);

        p1 = ro.add(rd.multiply(t1));
        p2 = ro.add(rd.multiply(t2));
        
        t = -1;
        if (t1 > 0 && withinReferenceCube(p1)) t = t1;
        else if (t2 > 0 && withinReferenceCube(p2)) t = t2;
        return t;
    }
    Vector3D getNormalAt (Vector3D intersectionPoint) override{
        double x, y, z, nx, ny, nz;
        Vector3D n;
        x = intersectionPoint.getX();
        y = intersectionPoint.getY();
        z = intersectionPoint.getZ();

        nx = 2 * A * x + D * y + E * z + G;
        ny = 2 * B * y + D * x + F * z + H;
        nz = 2 * C * z + E * x + F * y + I;

        n = Vector3D(nx, ny, nz);
        n.normalize();
        return n;
    }
};



class Floor : public Object {
public:
    Floor(double floorWidth, double tileWidth) {
        cout<<"Floor constructor"<<endl;
        cout<<"Floor width: "<<floorWidth<<endl;
        cout<<"Tile width: "<<tileWidth<<endl;
        reference_point = {- floorWidth / 2, - floorWidth / 2, 0};
        //print reference point
        cout<<"Reference point: "<<reference_point.getX()<<", "<<reference_point.getY()<<", "<<reference_point.getZ()<<endl;
        length = tileWidth;
    }

    void draw() override{
        int limit, i, j;
        bool isEven;
        double x, y;

        limit = -(int) (reference_point.getX());
        limit = limit  / length;

        glBegin(GL_QUADS);
        {
            for (i = -limit; i < limit; i++) {
                for (j = -limit; j < limit; j++) {
                    isEven = ((i + j) % 2 == 0);
                    
                    if (isEven) glColor3f(1.0, 1.0, 1.0);
                    else glColor3f(0.0, 0.0, 0.0);
                    
                    x = i * length;
                    y = j * length;
                    glVertex3f(x, y, 0);
                    glVertex3f(x + length, y, 0);
                    glVertex3f(x + length, y + length, 0);
                    glVertex3f(x, y + length, 0);
                }
            }
        }
        glEnd();
    }

    double intersect (Ray r, Color clr, int level , bool is_print = false) override{
        Vector3D n, ro, rd, intersectionPoint;
        double t;

        n = Vector3D(0, 0, 1);
        ro = r.getOrigin();
        rd = r.getDirection();

        t = (-1) * (n.dotProduct(ro) / n.dotProduct(rd));
        // check if intersection point lies within the floor
        intersectionPoint = ro.add(rd.multiply(t));

        if (intersectionPoint.getX() < reference_point.getX()) return -1;
        if (intersectionPoint.getX() > -reference_point.getX()) return -1;
        if (intersectionPoint.getY() < reference_point.getY()) return -1;
        if (intersectionPoint.getY() > -reference_point.getY()) return -1;
        return t;
    }

    Vector3D getNormalAt (Vector3D intersectionPoint) override{
        return Vector3D(0, 0, 1);
    }

    Color getColorAt (Vector3D intersectionPoint) override{
        int i,j;
        i = (reference_point.getX() + intersectionPoint.getX()) / length;
        j = (reference_point.getY() + intersectionPoint.getY()) / length;

        bool isEven = ((i + j) % 2 == 0);
                    
        if (isEven) return {1, 1, 1};
        return {0, 0, 0};
    }
};





#endif // RAYTRACING_1705051_HPP_INCLUDED
