#pragma once

#include <vector>
#include <glm/glm.hpp>



#include "LittleEngine/Graphics/Renderer.h"



namespace game
{
	//Simulates a 2D fluid in a grid using an Eulerian solver.
	class Fluid2D
	{

	public:

		Fluid2D() {};
		~Fluid2D() { Shutdown(); };

		bool Initialize(int width, int height, float h = 1.f, float viscosity = 0.001f);
		void Shutdown();
		void UpdateTimeStep(float dt);
		void Render(LittleEngine::Graphics::Renderer* renderer);
		void AddSource(int x, int y, float amount);
		void AddForce(int x, int y, float amountX, float amountY);

	private:

		void GetForcesAndSources();
		void UpdateFluidVelocityField();
		void UpdateFluidDensityField();


		bool IsBoundary(int i, int j);

		glm::vec2 TraceParticle(const glm::vec2& pos, float dt, const std::vector<glm::vec2>& velocityField);

#pragma region Velocity Field Methods

		glm::vec2 InterpolateVelocity(const glm::vec2& pos, const std::vector<glm::vec2>& velocityField);

		void AddVelocityForces();
		void AdvectVelocityField();
		void DiffuseVelocityField();
		void ProjectVelocityField();

#pragma endregion

#pragma region Density Field Methods

		float InterpolateDensity(const glm::vec2& pos, const std::vector<float>& densityField);

		void AddDensityForces();
		void AdvectDensityField();
		void DiffuseDensityField();
		void DissipateDensityField();

#pragma endregion

		int m_width = 0;
		int m_height = 0;

		float m_h = 1.f;						// grid cell size
		float m_dt = 0.1f;						// time step


		float m_viscosity = 0.001f;				// viscosity of the fluid


		std::vector<glm::vec2> m_velocity;		// velocity field
		std::vector<float> m_density;			// density field

		std::vector<glm::vec2> m_forces;
		std::vector<float> m_sources;		// density sources


		// intermediate fields
		std::vector<glm::vec2> w0;
		std::vector<glm::vec2> w1;
		std::vector<glm::vec2> w2;
		std::vector<glm::vec2> w3;

		std::vector<float> d0;
		std::vector<float> d1;
		std::vector<float> d2;
		std::vector<float> d3;


	};
}