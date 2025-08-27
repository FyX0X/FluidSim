#pragma once

#include <glm/glm.hpp>

namespace game
{


	struct GameData
	{
		glm::vec2 rectPos = { 1,1 };
		glm::vec2 pos2 = { 0,0 };
		float zoom = 50.f;

		

		// edge 1
		glm::vec2 A = { 0.f, 0.f };
		glm::vec2 B = { 1.f, 0.f };
		
		// edge 2
		glm::vec2 C = { 0.f, 1.f };
		glm::vec2 D = { 1.f, 1.f };


		LittleEngine::Graphics::Color color = LittleEngine::Graphics::Colors::Green;

	};

}