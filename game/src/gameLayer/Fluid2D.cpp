#include "Fluid2D.h"

#include <LittleEngine/Utils/logger.h>
#include <functional>

#include <Eigen/Sparse>

#include <iostream>

namespace game
{
	bool Fluid2D::Initialize(int width, int height, float h, float viscosity)
	{
		if (width <= 0 || height <= 0)
		{
			LittleEngine::Utils::Logger::Error("Fluid2D::Initialize: width and height must be greater than 0");
			return false;
		}

		if (h <= 0)
		{
			LittleEngine::Utils::Logger::Error("Fluid2D::Initialize: h must be greater than 0");
			return false;
		}


		m_h = h;

		m_width = width;
		m_height = height;

		m_viscosity = viscosity;

		m_forces.resize(m_width * m_height, glm::vec2(0.f, 0.f));
		m_sources.resize(m_width * m_height, 0.f);

		m_velocity.resize(m_width * m_height, glm::vec2(0.f, 0.f));
		w0.resize(m_width * m_height, glm::vec2(0.f, 0.f));
		w1.resize(m_width * m_height, glm::vec2(0.f, 0.f));
		w2.resize(m_width * m_height, glm::vec2(0.f, 0.f));
		w3.resize(m_width * m_height, glm::vec2(0.f, 0.f));

		m_density.resize(m_width * m_height, 0.f);
		d0.resize(m_width * m_height, 0.f);
		d1.resize(m_width * m_height, 0.f);
		d2.resize(m_width * m_height, 0.f);
		d3.resize(m_width * m_height, 0.f);


		std::cout << "sizes:\n";
		std::cout << m_velocity.size() << "\n";
		std::cout << w0.size() << "\n";
		std::cout << w1.size() << "\n";
		std::cout << m_forces.size() << "\n";
		std::cout << m_sources.size() << "\n";

		return true;
	}


	void Fluid2D::Shutdown()
	{
		m_h = 1.f;
		m_dt = 0.f;

		m_width = 0;
		m_height = 0;

		m_velocity.clear();
		m_density.clear();
	}

	void Fluid2D::UpdateTimeStep(float dt)
	{
		m_dt = dt;

		GetForcesAndSources();

		w0 = m_velocity;
		d0 = m_density;
		UpdateFluidVelocityField();
		UpdateFluidDensityField();


	}


	void Fluid2D::Render(LittleEngine::Graphics::Renderer* renderer)
	{

		// render density field as quads
		float scale = m_h * 0.5f;
		for (size_t j = 0; j < m_height; j++)
		{
			for (size_t i = 0; i < m_width; i++)
			{
				float density = m_density[j * m_width + i];
				if (density > 0.01f)
				{
					glm::vec4 color = glm::vec4(1.f, 1.f, 1.f, glm::clamp(density, 0.f, 1.f));
					renderer->DrawRect(glm::vec4{ (i + 0.5f) * m_h, (j + 0.5f) * m_h , scale, scale }, color);
				}
			}
		}


	}


	void Fluid2D::AddSource(int x, int y, float amount)
	{
		if (x < 0 || x >= m_width || y < 0 || y >= m_height)
			return;
		m_sources[y * m_width + x] += amount;
	}

	void Fluid2D::AddForce(int x, int y, float amountX, float amountY)
	{
		if (x < 0 || x >= m_width || y < 0 || y >= m_height)
			return;
		m_forces[y * m_width + x] += glm::vec2(amountX, amountY);
	}

	void Fluid2D::GetForcesAndSources()
	{
		m_forces.resize(m_width * m_height, glm::vec2(0.f, 0.f));
		m_sources.resize(m_width * m_height, 0.f);
		/// TODO: implement external forces and sources
	}


#pragma region Velocity Field Methods

	/**
	 * Calculates the fluid velocity field using the current velocity field and the forces acting on the fluid.
	 *
	 * add forces -> advect velocity field -> diffuse -> project to make it divergence free
	 *
	 */
	void Fluid2D::UpdateFluidVelocityField()
	{
		AddVelocityForces();
		AdvectVelocityField();
		DiffuseVelocityField();
		ProjectVelocityField();
	}

	void Fluid2D::AddVelocityForces()
	{
		for (size_t i = 0; i < m_velocity.size(); i++)
		{
			w1[i] = w0[i] + m_dt * m_forces[i];
		}
	}

	// pos is in world cooridinates ( 0 <= pos.x <= width*h, 0 <= pos.y <= height*h )
	glm::vec2 Fluid2D::InterpolateVelocity(const glm::vec2& pos, const std::vector<glm::vec2>& velocityField)
	{
		// bilinear interpolation to get velocity at position pos
		int x0 = static_cast<int>(pos.x / m_h);
		int y0 = static_cast<int>(pos.y / m_h);
		int x1 = x0 + 1;
		int y1 = y0 + 1;
		if (x0 < 0 || x1 >= m_width || y0 < 0 || y1 >= m_height)
			return glm::vec2(0.f);
		float sx = pos.x / m_h - x0;
		float sy = pos.y / m_h - y0;
		glm::vec2 v00 = velocityField[y0 * m_width + x0];
		glm::vec2 v10 = velocityField[y0 * m_width + x1];
		glm::vec2 v01 = velocityField[y1 * m_width + x0];
		glm::vec2 v11 = velocityField[y1 * m_width + x1];
		glm::vec2 v0 = (1 - sx) * v00 + sx * v10;
		glm::vec2 v1 = (1 - sx) * v01 + sx * v11;
		return (1 - sy) * v0 + sy * v1;
	}

	glm::vec2 Fluid2D::TraceParticle(const glm::vec2& pos, float dt, const std::vector<glm::vec2>& velocityField)
	{
		// RK2 integration to trace the particle backwards in time

		glm::vec2 k1 = InterpolateVelocity(pos, velocityField);
		glm::vec2 k2 = InterpolateVelocity(pos + dt * k1, velocityField);

		return pos + dt * 0.5f * (k1 + k2);
	}


	void Fluid2D::AdvectVelocityField()
	{
		// w2 = w1(p(X, - dt))		trace back the position of the fluid particle at time t-dt

		for (size_t j = 0; j < m_height; j++)
		{
			for (size_t i = 0; i < m_width; i++)
			{
				glm::vec2 pos = glm::vec2(i + 0.5f, j + 0.5f) * m_h;
				glm::vec2 tracedPos = TraceParticle(pos, -m_dt, w1);
				w2[j * m_width + i] = InterpolateVelocity(tracedPos, w1);
			}
		}
	}


	bool Fluid2D::IsBoundary(int i, int j)
	{
		return (i == 0 || i == m_width - 1 || j == 0 || j == m_height - 1);
	}


	void Fluid2D::DiffuseVelocityField()
	{
		// (I−nu * dt grad^2)w3(x) = w2(x)
		// solve using sparse methods (Ax = b).

		// TODO check if correct, waybe has to be done component wise.
		

		// Construct sparse matrix A and vectors bu, bv

		float coeff = m_viscosity * m_dt / (m_h * m_h);		
		std::vector<Eigen::Triplet<float>> tripletList;		// triplet is (row, col, value)
		Eigen::VectorXf bu(m_width * m_height);
		Eigen::VectorXf bv(m_width * m_height);
		for (size_t j = 0; j < m_height; j++)
		{
			for (size_t i = 0; i < m_width; i++)
			{
				int index = j * m_width + i;
				if (IsBoundary(i, j))
				{
					tripletList.emplace_back(index, index, 1); // boundary condition
					bu(index) = 0.f;
					bv(index) = 0.f;
				}

				tripletList.emplace_back(index, index, 1 + 4 * coeff); // center
				if (i > 1)				tripletList.emplace_back(index, index - 1, -coeff); // left
				if (i + 1 < m_width)	tripletList.emplace_back(index, index + 1, -coeff); // right
				if (j > 1)				tripletList.emplace_back(index, index - m_width, -coeff); // down
				if (j + 1 < m_height)	tripletList.emplace_back(index, index + m_width, -coeff); // up

				bu(index) = w2[index].x;
				bv(index) = w2[index].y;
			}
		}
		Eigen::SparseMatrix<float> A(m_width * m_height, m_width * m_height);
		A.insertFromTriplets(tripletList.begin(), tripletList.end());


		// solve Ax = b for each component

		Eigen::SimplicialLLT<Eigen::SparseMatrix<float>> solver;
		solver.compute(A);
		Eigen::VectorXf xu = solver.solve(bu);
		Eigen::VectorXf xv = solver.solve(bv);

		for (size_t j = 0; j < m_height; j++)
		{
			for (size_t i = 0; i < m_width; i++)
			{
				size_t index = j * m_width + i;
				w3[index] = glm::vec2(xu[index], xv[index]);
			}
		}


	}

	void Fluid2D::ProjectVelocityField()
	{
		// make the velocity field divergence free
		// solve Poisson equation using sparse methods (Ax = b).
		// w4(x) = w3(x)− grad q	(w4 is m_velocity, the final velocity field)
		// ∇^2q = ∇ • w3

		/*
		solve for q
		calculate ∇q
		calculate w4 = w3 - ∇q

		how to solve ∇²x = b:

		stencil for ∇²x:
			0  x  0
			x  ab x
			0  x  0
		/ beta

		beta = 1 + 4*a
		a = h^2 / (nu * dt)
		*/

		float alpha = (m_h * m_h) / (m_viscosity * m_dt);
		float beta = 1 + 4 * alpha;

		// Construct sparse matrix A and vector b
		std::vector<Eigen::Triplet<float>> tripletList;		// triplet is (row, col, value)
		Eigen::VectorXf b(m_width * m_height);

		for (size_t j = 0; j < m_height; j++)
		{
			for (size_t i = 0; i < m_width; i++)
			{
				size_t index = j * m_width + i;

				float neighborsSum = 0.f;

				if (i > 1)				{ tripletList.emplace_back(index, index - 1, -1.f);			neighborsSum += 1.f;} // left
				if (i + 1 < m_width)	{ tripletList.emplace_back(index, index + 1, -1.f);			neighborsSum += 1.f;} // right
				if (j > 1)				{ tripletList.emplace_back(index, index - m_width, -1.f);	neighborsSum += 1.f;} // down
				if (j + 1 < m_height)	{ tripletList.emplace_back(index, index + m_width, -1.f);	neighborsSum += 1.f; } // up
				tripletList.emplace_back(index, index, neighborsSum); // center

				float div = 0.f;
				// calculate divergence at (i, j) but skip boundaries
				if (i == 0 || i == m_width - 1 || j == 0 || j == m_height - 1)
				{
					b(index) = 0;
					continue;
				}
				float dudx = (w3[index+1].x - w3[index-1].y) / (2 * m_h);
				float dvdy = (w3[index + m_width].y - w3[index - m_width].y) / (2 * m_h);
				b(index) = -div;
			}
		}
		
		// solve Ax = b for q
		Eigen::SparseMatrix<float> A(m_width * m_height, m_width * m_height);
		A.insertFromTriplets(tripletList.begin(), tripletList.end());
		Eigen::SimplicialLLT<Eigen::SparseMatrix<float>> solver;
		solver.compute(A);
		Eigen::VectorXf q = solver.solve(b);
		// calculate ∇q and w4 = w3 - ∇q
		for (size_t j = 0; j < m_height; j++)
		{
			for (size_t i = 0; i < m_width; i++)
			{
				size_t index = j * m_width + i;
				if (IsBoundary(i, j))
				{
					m_velocity[index] = glm::vec2(0.f);
					continue;
				}
				float dqx = 0.f;
				float dqy = 0.f;
				if (i > 0 && i < m_width - 1)
					dqx = (q[index + 1] - q[index - 1]) / (2 * m_h);
				if (j > 0 && j < m_height - 1)
					dqy = (q[index + m_width] - q[index - m_width]) / (2 * m_h);
				m_velocity[index] = w3[index] - glm::vec2(dqx, dqy);
			}
		}


	}



#pragma endregion


#pragma region Density Field Methods


	void Fluid2D::UpdateFluidDensityField()
	{
		AddDensityForces();
		AdvectDensityField();
		DiffuseDensityField();
		DissipateDensityField();
	}

	void Fluid2D::AddDensityForces()
	{
		for (size_t i = 0; i < m_velocity.size(); i++)
		{
			w1[i] = w0[i] + m_dt * m_sources[i];
		}
	}

	float Fluid2D::InterpolateDensity(const glm::vec2& pos, const std::vector<float>& densityField)
	{
		// bilinear interpolation to get density at position pos
		int x0 = static_cast<int>(pos.x / m_h);
		int y0 = static_cast<int>(pos.y / m_h);
		int x1 = x0 + 1;
		int y1 = y0 + 1;
		if (x0 < 0 || x1 >= m_width || y0 < 0 || y1 >= m_height)
			return 0.f;
		float sx = pos.x / m_h - x0;
		float sy = pos.y / m_h - y0;
		float d00 = densityField[y0 * m_width + x0];
		float d10 = densityField[y0 * m_width + x1];
		float d01 = densityField[y1 * m_width + x0];
		float d11 = densityField[y1 * m_width + x1];
		float d0 = (1 - sx) * d00 + sx * d10;
		float d1 = (1 - sx) * d01 + sx * d11;
		return (1 - sy) * d0 + sy * d1;
	}

	void Fluid2D::AdvectDensityField()
	{
		// d2 = d1(p(X, - dt))		trace back the position of the fluid particle at time t-dt
		for (size_t j = 0; j < m_height; j++)
		{
			for (size_t i = 0; i < m_width; i++)
			{
				glm::vec2 pos = glm::vec2(i + 0.5f, j + 0.5f) * m_h;
				glm::vec2 tracedPos = TraceParticle(pos, -m_dt, m_velocity);
				d2[j * m_width + i] = InterpolateDensity(tracedPos, d1);
			}
		}
	}

	void Fluid2D::DiffuseDensityField()
	{
		// (I−nu * dt grad^2)d3(x) = d2(x)
		// solve using sparse methods (Ax = b).
		float coeff = m_viscosity * m_dt / (m_h * m_h);
		std::vector<Eigen::Triplet<float>> tripletList;		// triplet is (row, col, value)
		Eigen::VectorXf b(m_width * m_height);
		for (size_t j = 0; j < m_height; j++)
		{
			for (size_t i = 0; i < m_width; i++)
			{
				int index = j * m_width + i;
				if (IsBoundary(i, j))
				{
					tripletList.emplace_back(index, index, 1); // boundary condition
					b(index) = 0;
				}
				tripletList.emplace_back(index, index, 1 + 4 * coeff); // center
				if (i > 1)				tripletList.emplace_back(index, index - 1, -coeff); // left
				if (i + 1 < m_width)	tripletList.emplace_back(index, index + 1, -coeff); // right
				if (j > 1)				tripletList.emplace_back(index, index - m_width, -coeff); // down
				if (j + 1 < m_height)	tripletList.emplace_back(index, index + m_width, -coeff); // up
				b(index) = d2[index];
			}
		}
		Eigen::SparseMatrix<float> A(m_width * m_height, m_width * m_height);
		A.insertFromTriplets(tripletList.begin(), tripletList.end());

		// solve Ax = b
		Eigen::SimplicialLLT<Eigen::SparseMatrix<float>> solver;
		solver.compute(A);
		Eigen::VectorXf x = solver.solve(b);
		for (size_t j = 0; j < m_height; j++)
		{
			for (size_t i = 0; i < m_width; i++)
			{
				size_t index = j * m_width + i;
				d3[index] = x[index];
			}
		}
	}


	void Fluid2D::DissipateDensityField()
	{
		// (1+∆tα)d4(x,t +∆t) = d3(x,t)

		float alpha = 0.01f; // dissipation rate
		float coeff = 1.f / (1.f + m_dt * alpha);
		for (size_t i = 0; i < m_density.size(); i++)
		{
			m_density[i] = coeff * d3[i];
		}
	}

#pragma endregion



	/**
	 * Fourth-order Runge-Kutta method for solving ordinary differential equations (ODEs).
	 * calculates y(t + dt) given y(t) and the derivative function f.
	 * f(t) = dy/dt
	 */
	float RungeKutta4(float y, float dt, std::function<float(float)> f)
	{

		float k1 = f(y);
		float k2 = f(y + 0.5f * dt * k1);
		float k3 = f(y + 0.5f * dt * k2);
		float k4 = f(y + dt * k3);
		return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0f;

	}

}
