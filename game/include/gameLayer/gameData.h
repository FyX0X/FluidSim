#pragma once

#include <glm/glm.hpp>

#include "Fluid2D.h"

namespace game
{


	struct GameData
	{
		glm::vec2 cameraPos = { 0.f, 0.f };
		Fluid2D fluid;

	};

}