# ASE Calculator Builders

One of the key parts of the [`ASENEBCalculator`](@ref) kinetic calculator is the idea of an energy/force calculator *builder*. These are Julia functors that instantiate the Python `Calculator` classes responsible for determining the energy of a given atomic system in ASE. These builders are required because it is often necessary to construct multiple ASE `Calculator`s for different `Atoms` objects, but there is no simple way to do this automatically within Kinetica's calculation workflows. Builders therefore abstract this construction away from the user, presenting a handful of parameters with sensible defaults for a given calculator.

Kinetica.jl comes with three builders, but they're very easy to implement - not only for `Calculator`s that exist within ASE, but for others that come from external packages and themselves have interfaces to ASE. For example, while ASE doesn't internally support any MLIPs for energy and force evaluation, many MLIPs come with Python frontends that can hook into ASE. Kinetica can take advantage of these frontends to also hook such MLIPs into its rate calculation workflow.

As an example, this tutorial will implement a calculator builder for [NeuralNEB](https://arxiv.org/abs/2207.09971), a modern MLIP based on the PaiNN architecture. NeuralNEB's [Python frontend](https://gitlab.com/matschreiner/neuralneb) comes with an ASE calculator under `neuralneb.utils.MLCalculator`, which requires a PaiNN model (which can be constructed through `neuralneb.painn.PaiNN`) that has been loaded with trained parameters as input.

At their core, most calculator builders have 3 components:

* A struct containing parameters.
* An outer constructor providing a way of loading Python objects into the struct once at its construction. This is not technically necessary, but speeds up the creation of ASE energy/force calculators on-the-fly since it avoids excessive reimporting of Python modules.
* A functor method with arguments `(dir::String, mult::Int, chg::Int, kwargs...)`, where `dir` is the directory a given calculation should take place in, `mult` is the spin multiplicity of a supplied atomic system, `chg` is its charge and `kwargs` are keyword arguments to pass to the ASE calculator class upon instantiation. Some or all of these arguments may not be used depending on the calculator, but Kinetica will always provide them.

For our exemplary NeuralNEB calculator, its struct could look something like this (note that PythonCall needs to be accessible):

```julia
mutable struct NeuralNEBBuilder
    calc_class::Py
    painn_class::Py
    device::String
    statedict
end
```

Here, `calc_class` will hold the uninstantiated ASE calculator class, `painn_class` will do the same for the uninstantiated PaiNN class, `device` will be the PyTorch device to run the calculations on (CPU or GPU) and `statedict` will hold the parameters for the trained model. Rather than users having to pass this all manually and import these classes themselves though, we can implement an outer constructor:

```julia
function NeuralNEBBuilder(modelpath::String, device::String)
    torch = pyimport("torch")
    nn_utils = pyimport("neuralneb.utils")
    nn_painn = pyimport("neuralneb.painn")
    statedict = torch.load(modelpath)
    return NeuralNEBBuilder(nn_utils.MLCalculator, nn_painn.PaiNN, device, statedict)
end
```

Now all that users of our builder need to supply are the path to the model parameters (`modelpath`) and the device. The constructor does the heavy lifting of setting up the Python side of things.

Finally, this struct needs a method that can be called to instantiate the `calc_class`, such that it can be provided as a calculator directly to ASE:

```julia
function (builder::NeuralNEBBuilder)(dir::String, mult::Int, chg::Int, kwargs...)
    model = builder.painn_class(3, 256, 5)
    model.load_state_dict(builder.statedict)
    model.eval()
    return builder.calc_class(model, device=builder.device)
end
```

This functor creates a fresh PaiNN model from the parameters in the struct, sets it to evaluation mode and instantiates the `MLCalculator` class on the requested device. Notably, none of the arguments of this functor are actually used. This is fine, as NeuralNEB doesn't write any files, and doesn't take any external multiplicity or charge information. Many other calculators will do so however, so these arguments are provided as a catch-all to ensure that all builders have the information required to correctly instantiate their calculators.

After that, the builder is complete! It can now be used like any other energy/force calculator currently in Kinetica by constructing the builder and passing it to the [`ASENEBCalculator`](@ref):

```julia
builder = NeuralNEBBuilder(model="/path/to/model.sd", device="cuda")
calc = ASENEBCalculator(builder, "./calc")

# solvemethod = VariableODESolve(pars, conditions, calc)
# res = solve_network(...)
```