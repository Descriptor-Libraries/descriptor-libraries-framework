Substructure Search
====================

The substructure search page allows the user to search the database using a chemical substructure.
There are two methods users can use to input their substructure: (1) Drawing the molecule using the molecular sketcher
and (2) Inputting a SMILES or SMARTS string into the search box.

General Approach
----------------
This page utilizes state variables to define search parameters, such as the search target (smiles or smarts string), search limit, and offset. 
A search is triggered whenever the `searchToggle` state variable changes. 
Currently, the ``searchToggle`` is called in the functions ``newSearch`` and ``loadMore``, both of these functions update state variables for search. 
These are the functions that should be called if a new search is started or more results are loaded.
When a new search is performed, search variables like limit and offset are set to default values. 
If more results are loaded, the state variables are incremented accordingly, and the relevant parts of the page are updated.

The pages is updated through an effect hook that calls a function called ``loadImages``

.. code-block:: javascript

    // initial load of data
    // and load when search changes. 
    useEffect( ( ) => {
        const controller = new AbortController();
        const signal = controller.signal;

        loadImages(signal);
        
        return () => {
          controller.abort();
        }
      },
        // eslint-disable-next-line react-hooks/exhaustive-deps
        [ searchToggle ] 
    );


The ``loadImages`` makes two API calls.
First, the search results are retrieved from the 
database using the state varaibles as inputs to the `search` API endpoint. The SMILES strings from these results are then used
to retrieve images from the `depict` endpoint. The images are then populated into a component called ``DynamicGric.``

Updating the page is handled using the hook shown on this page and happens any time the variable ``searchToggle`` changes.
One variable has been chosen to represent a search change. Using all of the search variables as dependencies in ``useEffect`` would 
result in multiple page loads if more than one search parameters changes. 
``searchToggle`` is called in the ``newSearch`` and ``loadMore`` functions.
The ``loadMore`` function increments search parameters limit and offset in order to load more on the page.
The ``newSearch`` function updates the search string and sets search parameters to default values. 

Room for Improvement
--------------------
This code should probably at a minimum be refactored to better separate component and search responsibility.
In the future, we should consider adding more advanced state management.
It's functional for now though!


