#include <algorithm>
#include <iostream>
#include <vector>
#include <random>
#include <limits>
#include <SFML/graphics.hpp>

#include "HSL.hpp"


#define SQ(x) (x * x)

enum class ScalarField
{
    velocity,
    pressure
};

auto display = ScalarField::pressure;

float rho = 1.0f;
float vis = 0.5f;
float gx = 0.0f;
float gy = -10.0f;

float domainX = 2.0f;
float domainY = 2.0f;

unsigned int nx = 81u;
unsigned int ny = 81u;

float dx = domainX / (nx - 1);
float dy = domainY / (ny - 1);
float dt = 0.0001;

unsigned int windowWidth = 800;
unsigned int windowHeight = 800;

float pixelX = static_cast<float>(windowWidth) / static_cast<float>(nx - 1);
float pixelY = static_cast<float>(windowHeight) / static_cast<float>(ny - 1);

std::vector<std::vector<float>> sourceField;

std::vector<std::vector<float>> pressureField;
std::vector<std::vector<float>> nextPressureField;

std::vector<std::vector<sf::Vector2f>> velocityField;
std::vector<std::vector<sf::Vector2f>> nextVelocityField;

unsigned int fieldLineNumber = 20;
std::vector<std::pair<sf::Vector2f, unsigned int>> velocityFieldLines;

int counter = 0;


void Init ()
{
    sourceField.resize(nx);
    pressureField.resize(nx);
    velocityField.resize(nx);

    for (int i = 0; i < nx; ++i)
    {
        sourceField[i].resize(ny, 0.0f);
        pressureField[i].resize(ny, 1.0f);
        velocityField[i].resize(ny, sf::Vector2f(2.0f, 0.0f));
    }

    nextPressureField = pressureField;
    nextVelocityField = velocityField;
}


void BoundaryCondition (float time)
{
    // bottom and top
    for (int i = 0; i < nx; ++i)
    {
        pressureField[i][0] = pressureField[i][1];
        pressureField[i][ny - 1] = 1.0f;

        velocityField[i][0].x = 0.0f;
        velocityField[i][0].y = 0.0f;
        // velocityField[i][ny - 1].x = cos(25.0f * time);
        velocityField[i][ny - 1].x = 1.0;
        velocityField[i][ny - 1].y = 0.0f;
    }

    // left and right
    for (int j = 0; j < ny; ++j)
    {
        pressureField[0][j] = pressureField[1][j];
        pressureField[nx - 1][j] = pressureField[nx - 2][j];

        velocityField[0][j].x = 0.0f;
        velocityField[0][j].y = 0.0f;
        velocityField[nx - 1][j].x = 0.0f;
        velocityField[nx - 1][j].y = 0.0f;
    }
}


void UpdatePressureFieldSource ()
{
    for (int i = 1; i < (nx - 1); ++i)
    {
        for (int j = 1; j < (ny - 1); ++j)
        {
            auto dudx =
                (velocityField[i + 1][j].x - velocityField[i - 1][j].x)
                / (2.0f * dx);
            auto dudy =
                (velocityField[i][j + 1].x - velocityField[i][j - 1].x)
                / (2.0f * dy);
            auto dvdx =
                (velocityField[i + 1][j].y - velocityField[i - 1][j].y)
                / (2.0f * dx);
            auto dvdy =
                (velocityField[i][j + 1].y - velocityField[i][j - 1].y)
                / (2.0f * dy);
            auto delU = dudx + dvdy;
            sourceField[i][j] =
                rho * (delU / dt + SQ(dudx) + 2 * dudy * dvdx + SQ(dvdy));
        }
    }
}


void UpdatePressureField (int iterations)
{
    for (int it = 0; it < iterations; ++it)
    {
        for (int i = 1; i < (nx - 1); ++i)
        {
            for (int j = 1; j < (ny - 1); ++j)
            {
                auto px = pressureField[i + 1][j] + pressureField[i - 1][j];
                auto py = pressureField[i][j + 1] + pressureField[i][j - 1];
                auto b = sourceField[i][j];

                nextPressureField[i][j] =
                    (SQ(dy) * px + SQ(dx) * py - b * SQ(dx) * SQ(dy)) /
                    (2 * (SQ(dx) + SQ(dy)));
            }
        }

        pressureField = nextPressureField;
        BoundaryCondition(static_cast<float>(counter) * dt);
    }
}


void UpdateVelocityField ()
{
    for (int i = 1; i < (nx - 1); ++i)
    {
        for (int j = 1; j < (ny - 1); ++j)
        {
            auto u = velocityField[i][j].x;
            auto v = velocityField[i][j].y;

            auto dudx =
                (velocityField[i + 1][j].x - velocityField[i - 1][j].x)
                / (2.0f * dx);
            auto dudy =
                (velocityField[i][j + 1].x - velocityField[i][j - 1].x)
                / (2.0f * dy);
            auto dvdx =
                (velocityField[i + 1][j].y - velocityField[i - 1][j].y)
                / (2.0f * dx);
            auto dvdy =
                (velocityField[i][j + 1].y - velocityField[i][j - 1].y)
                / (2.0f * dy);

            auto dudxdx =
                (velocityField[i + 1][j].x - 2.0f * velocityField[i][j].x +
                 velocityField[i - 1][j].x)
                / SQ(dx);
            auto dudydy =
                (velocityField[i][j + 1].x - 2.0f * velocityField[i][j].x +
                 velocityField[i][j - 1].x)
                / SQ(dy);
            auto dvdxdx =
                (velocityField[i + 1][j].y - 2.0f * velocityField[i][j].y +
                 velocityField[i - 1][j].y)
                / SQ(dx);
            auto dvdydy =
                (velocityField[i][j + 1].y - 2.0f * velocityField[i][j].y +
                 velocityField[i][j - 1].y)
                / SQ(dy);

            auto dpdx = (pressureField[i + 1][j] - pressureField[i - 1][j])
                / (2.0f * dx);
            auto dpdy = (pressureField[i][j + 1] - pressureField[i][j - 1])
                        / (2.0f * dy);

            nextVelocityField[i][j].x = u + dt * (vis * (dudxdx + dudydy) + gx - (dpdx / rho) - (u * dudx) - (v * dudy));
            nextVelocityField[i][j].y = v + dt * (vis * (dvdxdx + dvdydy) + gy - (dpdy / rho) - (u * dvdx) - (v * dvdy));
        }
    }

    velocityField = nextVelocityField;
    BoundaryCondition(static_cast<float>(counter) * dt);
}


sf::Vector2f GetVelocityAt (float x, float y)
{
    auto i = std::max(0u, std::min(static_cast<unsigned int>(x / dx), nx - 2));
    auto j = std::max(0u, std::min(static_cast<unsigned int>(y / dy), ny - 2));

    auto px = (x - i * dx) / dx;
    auto py = (y - j * dy) / dy;

    auto v00 = velocityField[i    ][j   ];
    auto v10 = velocityField[i + 1][j   ];
    auto v11 = velocityField[i + 1][j + 1];
    auto v01 = velocityField[i    ][j + 1];

    auto v0 = (1.0f - px) * v00 + px * v01;
    auto v1 = (1.0f - px) * v10 + px * v11;

    auto velocity = (1.0f - py) * v0 + py * v1;

    return velocity;
}


float GetMaxPressure ()
{
    float max = -100000.0f;

    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            max = std::max(max, pressureField[i][j]);
        }
    }

    return max;
}


float GetMinPressure ()
{
    float min = 100000.0f;

    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            min = std::min(min, pressureField[i][j]);
        }
    }

    return min;
}


float GetMaxSpeed ()
{
    float max = -1.0f;

    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            auto velocity = velocityField[i][j];
            max = std::max(max, sqrtf(SQ(velocity.x) + SQ(velocity.y)));
        }
    }

    return max;
}


sf::Vector2f GetMeshCoordinate(int i, int j)
{
    return sf::Vector2f(i * pixelX, (ny - j - 1) * pixelY);
}


sf::Vector2f GetDomainCoordinate(float x, float y)
{
    return sf::Vector2f(static_cast<float>(windowWidth) * x / domainX,
                        static_cast<float>(windowHeight) * (1.0f - y / domainY));
}


sf::Color CreateColor (float value)
{
    value = std::min(value, 1.0f);
    auto h = (1.0f - value) * 240.0f;
    auto hsl = HSL(static_cast<int>(h), 100, 50);
    return hsl.TurnToRGB();
}


sf::Color CreatePressureColor (int i, int j, float minPressure, float maxPressure)
{
    auto pressure = pressureField[i][j];
    auto range = (pressure - minPressure) / (maxPressure - minPressure);
    return CreateColor(range);
}


sf::Vertex CreatePressureVertex (int i, int j, float minPressure, float maxPressure)
{
    return sf::Vertex(GetMeshCoordinate(i, j), CreatePressureColor(i, j, minPressure, maxPressure));
}


sf::VertexArray CreatePressureCell (int i, int j, float minPressure, float maxPressure)
{
    sf::VertexArray va(sf::Quads);

    va.append(CreatePressureVertex(i,     j,     minPressure, maxPressure));
    va.append(CreatePressureVertex(i + 1, j,     minPressure, maxPressure));
    va.append(CreatePressureVertex(i + 1, j + 1, minPressure, maxPressure));
    va.append(CreatePressureVertex(i,     j + 1, minPressure, maxPressure));

    return va;
}


sf::Color CreateSpeedColor (int i, int j, float maxSpeed)
{
    auto velocity = velocityField[i][j] / maxSpeed;
    auto speed = sqrtf(SQ(velocity.x) + SQ(velocity.y));
    return CreateColor(speed);
}


sf::Vertex CreateSpeedVertex (int i, int j, float maxSpeed)
{
    return sf::Vertex(GetMeshCoordinate(i, j), CreateSpeedColor(i, j, maxSpeed));
}


sf::VertexArray CreateSpeedCell (int i, int j, float maxSpeed)
{
    sf::VertexArray va(sf::Quads);

    va.append(CreateSpeedVertex(i,     j,     maxSpeed));
    va.append(CreateSpeedVertex(i + 1, j,     maxSpeed));
    va.append(CreateSpeedVertex(i + 1, j + 1, maxSpeed));
    va.append(CreateSpeedVertex(i,     j + 1, maxSpeed));

    return va;
}


sf::VertexArray DrawFieldLine (sf::Vector2f position, unsigned int n)
{
    auto va = sf::VertexArray(sf::LineStrip);

    for (int it = 0; it < n; ++it)
    {
        va.append(sf::Vertex(GetDomainCoordinate(position.x, position.y), sf::Color::Black));
        auto velocity = GetVelocityAt(position.x, position.y);
        position += 0.001f * velocity;
    }

    return va;
}


void Draw (sf::RenderWindow& window)
{
    auto minPressure = GetMinPressure();
    auto maxPressure = GetMaxPressure();
    auto maxSpeed = GetMaxSpeed();

    auto va2 = sf::VertexArray(sf::Lines);

    for (int i = 0; i < (nx - 1); ++i)
    {
        for (int j = 0; j < (ny - 1); ++j)
        {
            sf::VertexArray va;
            switch (display)
            {
                case ScalarField::pressure:
                    va = CreatePressureCell(i, j, minPressure, maxPressure);
                    break;
                case ScalarField::velocity:
                    va = CreateSpeedCell(i, j, maxSpeed);
                    break;
            }
            window.draw(va);

            if (i % 4 == 2 && j % 4 == 2)
            {

                auto position = GetMeshCoordinate(i, j);
                auto velocity = GetVelocityAt(i * dx, j * dy);

                auto speed = sqrtf(SQ(velocity.x) + SQ(velocity.y));
                auto scaled = 10 * log(speed + 1.2f);
                velocity *= 10.0f * scaled / speed;
                velocity.y *= -1.0f;
                va2.append(sf::Vertex(position, sf::Color::Black));
                va2.append(sf::Vertex(position + velocity, sf::Color::Black));

                auto circle = sf::CircleShape(5.0f, 10);
                circle.move(position);
                circle.move(-5.0f, -5.0f);
                circle.setFillColor(sf::Color::Black);
                window.draw(circle);
            }
        }
    }

    window.draw(va2);
}


int main ()
{
    Init();
    BoundaryCondition(static_cast<float>(counter) * dt);

    sf::ContextSettings settings;
    settings.antialiasingLevel = 2;

    sf::RenderWindow window(sf::VideoMode(windowWidth, windowHeight),
                            std::string("Navier-Stokes Equations."),
                            sf::Style::Default,
                            settings);

    while (window.isOpen())
    {
        auto event = sf::Event();
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
            {
                window.close();
            }
            if (event.type == sf::Event::KeyPressed)
            {
                if (event.key.code == sf::Keyboard::Key::Escape)
                {
                    window.close();
                }
                if (event.key.code == sf::Keyboard::Key::D)
                {
                    display = static_cast<ScalarField>(1 - static_cast<int>(display));
                }
                if (event.key.code == sf::Keyboard::Key::S)
                {
                    dt = 0.0001;
                }
                if (event.key.code == sf::Keyboard::Key::R)
                {
                    dt /= 1.5f;
                }
                if (event.key.code == sf::Keyboard::Key::T)
                {
                    dt *= 1.5f;
                }
            }
        }

        switch (display)
        {
            case ScalarField::velocity:
                window.setTitle(std::string("Navier-Stokes Equations. Display Velocity. t = ")
                                + std::to_string(static_cast<float>(counter) * dt));
                break;
            case ScalarField::pressure:
                window.setTitle(std::string("Navier-Stokes Equations. Display Pressure. t = ")
                                + std::to_string(static_cast<float>(counter) * dt));
                break;
        }

        UpdatePressureFieldSource();
        UpdatePressureField(100);
        UpdateVelocityField();

        ++counter;

        Draw(window);

//        window.setFramerateLimit(300);
        window.display();
    }

    return 0;
}